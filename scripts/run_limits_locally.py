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

args = parser.parse_args()
index_list = range(0, 2820)


def command(index):
    os.system("./{} {}".format(args.task, str(index)))


p = Pool(args.ncores)
p.map(command, index_list)