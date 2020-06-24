#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import sys
from multiprocessing import Pool

def execute(cmd):
    try:
        os.system(cmd)
    except:
        print "[WARNING] Command failed:",cmd


datacard_pattern = sys.argv[1]
freezeargs = sys.argv[2]
wsname = sys.argv[3]
ncores = 10
if len(sys.argv) > 4:
    ncores = int(sys.argv[4])

datacards = [d for d in glob.glob(datacard_pattern) if "/cmb/" not in d]

cmds = ['card=DATACARD; basedir=$(dirname $(dirname ${card})); category=$(basename $(dirname ${card})); echo $card; echo ${card/combined.txt.cmb/WSNAME}; echo ${card/combined.txt.cmb/prefit_shapes.root}; PostFitShapesFromWorkspace -w ${card/combined.txt.cmb/WSNAME}  -o ${card/combined.txt.cmb/prefit_shapes.root} -d ${basedir}/restore_binning/${category}/${category}.txt FREEZEARGS'.replace("DATACARD",d).replace("FREEZEARGS",freezeargs).replace("WSNAME",wsname) for d in datacards]

p = Pool(ncores)
p.map(execute, cmds)
