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
print sys.argv
ws_pattern = sys.argv[1]
freezeargs = sys.argv[2]
ncores = 10
if len(sys.argv) > 3:
    ncores = int(sys.argv[3])

ws = [w for w in glob.glob(ws_pattern) if "/htt_" in w]

cmds = ['ws=WORKSPACE; wsname=$(basename ${ws}); basedir=$(dirname $(dirname ${ws})); category=$(basename $(dirname ${ws})); PreFitSignalShapes -w ${ws}  -o ${ws/${wsname}/prefit_signal_shapes.root} -d ${basedir}/restore_binning/${category}/${category}.txt -c ${category} --freeze FREEZEARGS'.replace("WORKSPACE",w).replace("FREEZEARGS",freezeargs) for w in ws]

p = Pool(ncores)
p.map(execute, cmds)
