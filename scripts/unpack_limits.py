#!/usr/bin/env python

import tarfile
import sys
import os
import glob

from multiprocessing import Pool


def extract(info):

    if tarfile.is_tarfile(info[0]):
        f = tarfile.open(info[0],"r")
        f.extractall(path=info[1])
        return 0
    else:
        return 1
    

pattern = sys.argv[1]
outfolder = sys.argv[2]

info_list = [(i, outfolder) for i in glob.glob(pattern)]

p = Pool(10)
returncodes = p.map(extract,info_list)
print("Sum of returncodes",sum(returncodes))
