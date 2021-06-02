import argparse
import os
import shutil
import tarfile


def parse_arguments():
    parser = argparse.ArgumentParser(description="Run combine fits remotely")
    parser.add_argument("--combine-script",
                        type=str,
                        required=True,
                        help="path to the combine script with the fits")
    parser.add_argument(
        "--workspace",
        type=str,
        required=True,
        help="path to the workspace used for the fit",
    )
    parser.add_argument(
        "--workdir",
        type=str,
        required=True,
        help="path to the location of the gc_workdir",
    )
    parser.add_argument(
        "--se-path",
        type=str,
        required=True,
        help="path where to store the samples on nrg",
    )
    parser.add_argument("--tag",
                        help="Tag for the gc task",
                        default="combine_fits")

    return parser.parse_args()


def filter_MSSMCode(tarinfo):
    if "CombineHarvester/MSSMvsSMRun2Legacy" in (tarinfo.name):
        return None
    return tarinfo


def parse_nrg_path(path):
    if path.startswith("/storage/gridka-nrg"):
        se_path = path.replace("/storage/gridka-nrg/",
                               "root://cmsxrootd-kit.gridka.de//store/user/")
    elif path.startswith(
            "srm://cmssrm-kit.gridka.de:8443/srm/managerv2?SFN=/pnfs/gridka.de/cms/disk-only/store/user/"
    ):
        se_path = path.replace(
            "srm://cmssrm-kit.gridka.de:8443/srm/managerv2?SFN=/pnfs/gridka.de/cms/disk-only/store/user/",
            "root://cmsxrootd-kit.gridka.de//store/user/")
    elif path.startswith("root://cmsxrootd-kit.gridka.de//store/user/"):
        se_path = path
    else:
        print("se-path unknown, aborting")
        raise Exception
    return se_path


def read_njobs(config):
    njobs = 0
    return njobs


def upload_tarball(se_path, extra_files, outputpath):
    print("building tarball...")
    outputfile = "cmssw.tar.gz"
    tar = tarfile.open(outputfile, "w:gz", dereference=True)
    tar.add(os.environ.get("CMSSW_BASE"),
            recursive=True,
            filter=filter_MSSMCode,
            arcname=os.environ.get("CMSSW_VERSION"))
    for file in extra_files:
        tar.add(file, os.path.basename(file))
    tar.close()
    print("finished building tarball...")
    print("upload tarball...")
    cmd = "xrdcp -fp {outputfile} {tarballpath}/{TARBALLNAME}".format(
        outputfile=outputfile, tarballpath=se_path, TARBALLNAME=outputfile)
    print(cmd)
    os.system(cmd)
    print("finished uploading tarball...")
    print("Create output folder")
    cmd = "xrdfs root://cmsxrootd-kit.gridka.de/ mkdir {outputpath}".format(
        outputpath=outputpath.replace("root://cmsxrootd-kit.gridka.de/", ""))
    os.system(cmd)
    return "{tarballpath}/{TARBALLNAME}".format(tarballpath=se_path,
                                                TARBALLNAME=outputfile)


def modify_combine_script(workspace, script, workdir):
    print("adapting combine config ...")
    modified_config = open("{}/combine_fit.sh".format(workdir), "w+")
    njobs = 0
    with open(script, "r") as file:
        for line in file:
            if " -eq " in line:
                njobs += 1
            if line.startswith("cd") or line.startswith(
                    "source") or line.startswith("eval"):
                continue
            elif workspace in line:
                modified_config.write(line.replace(workspace, "ws.root"))
            else:
                modified_config.write(line)
    modified_config.close()
    os.chmod("{}/combine_fit.sh".format(workdir), 0o777)
    print("Submitting {} jobs".format(njobs))
    return njobs


def write_gc(script, workspace, workdir, tag, se_path):
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    workspace = os.path.abspath(workspace)
    script = os.path.abspath(script)
    workdir = os.path.abspath(workdir)
    se_path = parse_nrg_path(se_path)
    print("Using se path: {}".format(se_path))

    configfilepath = "{WORKDIR}/{TAG}.conf".format(WORKDIR=workdir, TAG=tag)
    outputpath = '{se_path}/output/'.format(se_path=se_path)
    # copy config, worksapce and fitiing script to the correct location for easy zipping
    njobs = modify_combine_script(workspace, script, workdir)
    shutil.copy2("scripts/default_gc_etp7.conf", configfilepath)
    shutil.copy2("scripts/run_combine_remote.sh",
                 "{}/run_combine_remote.sh".format(workdir))
    shutil.copy2(workspace, "{}/ws.root".format(workdir))

    extra_files = [
        configfilepath, "{}/combine_fit.sh".format(workdir),
        "{}/ws.root".format(workdir)
    ]
    tarballpath = upload_tarball(se_path, extra_files, outputpath)
    configfile = open(configfilepath, "a+")

    configfile.write("\n[constants]\n")
    configfile.write(
        'TARBALL_PATH = {tarballpath}\n'.format(tarballpath=tarballpath))
    configfile.write(
        'OUTPUT_PATH = {outputpath}\n'.format(outputpath=outputpath))
    configfile.write('SCRAM_ARCH = {}\n'.format(os.environ.get("SCRAM_ARCH")))
    configfile.write('CMSSW_VERSION = {}\n'.format(
        os.environ.get("CMSSW_VERSION")))
    configfile.write("[parameters]\n")
    configfile.write('JOB_ID = range(0,{njobs})\n'.format(njobs=njobs))
    configfile.close()
    return "go {config} -G -m 0".format(config=configfilepath)


def main(args):
    # check if CMSSW is there:
    if "SCRAM_ARCH" in os.environ and "CMSSW_BASE" in os.environ and "CMSSW_VERSION" in os.environ:
        print(write_gc(args.combine_script, args.workspace, args.workdir, args.tag,
                        args.se_path))
    else:
        print("No CMSSW found ... Exiting ")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
