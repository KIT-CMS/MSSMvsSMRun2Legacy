import os
import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser(
    description="Run Model dependent limits locally")

parser.add_argument("--taskname",
                    type=str,
                    required=True,
                    help="name of the condor task script")
parser.add_argument('--cores', default=20, help="number of cores to be used")
parser.add_argument('--njobs', default=2820, help="number of jobs that are processed")

args = parser.parse_args()
index_list = range(0, int(args.njobs))


def command(index):
    os.system("./{} {}".format(args.taskname, str(index)))


p = Pool(int(args.cores))
p.map(command, index_list)