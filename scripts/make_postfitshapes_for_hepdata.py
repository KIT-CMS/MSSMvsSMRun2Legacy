#! /usr/bin/env python

from argparse import ArgumentParser
from multiprocessing import Pool
import yaml
import os
import subprocess
import shlex

def execute_command(command):
  print("Processing:",command)
  returncode = subprocess.check_call(shlex.split(command), stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
  return returncode

parser = ArgumentParser(description="Script to create postfit shapes for HEP Data using PostFitShapesForHEPData")
parser.add_argument("--analysis-configuration", required=True, help="Path to a .yaml file containing all information for an analysis.")

args = parser.parse_args()

command_template = ("PostFitShapesForHEPData"
  " -w {WORKSPACE}"
  " -d {RESTORE_BINNING_DATACARD}"
  " -c {CATEGORY}"
  " -f {FITFILE} -F {FITNAME}"
  " {POI_CONFIG}"
  " -m {MASSNAME} {MASSES}"
)

commands = []

analysis_configuration = yaml.load(open(args.analysis_configuration, "r"))

pois = " ".join(["-P {key}:{value}".format(key=key,value=value) for key, value in analysis_configuration["POI"].items()])
masses = " ".join(["-M {m}".format(m=m) for m in analysis_configuration["masses"]])

for era in analysis_configuration["eras"]:
  for fs in analysis_configuration["final_states"]:
    workspace = analysis_configuration["{fs}_{era}_workspace".format(fs=fs, era=era)]
    for c in analysis_configuration["{fs}_categories".format(fs=fs)]:
      cname = "_".join(["htt",fs,str(c),str(era)])
      restore_datacard = os.path.join(analysis_configuration["restore_directory"], cname+".txt")
      commands.append(command_template.format(
          WORKSPACE=workspace,
          RESTORE_BINNING_DATACARD=restore_datacard,
          CATEGORY=cname,
          FITFILE=analysis_configuration["fitfile"],
          FITNAME=analysis_configuration["fitname"],
          POI_CONFIG=pois,
          MASSNAME=analysis_configuration["massname"],
          MASSES=masses
        )
      )

p = Pool(10)
returncodes = p.map(execute_command, commands)
print("Sum of returncodes:",sum(returncodes))
