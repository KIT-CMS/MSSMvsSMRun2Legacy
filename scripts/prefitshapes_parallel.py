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
ncores = 10
if len(sys.argv) > 2:
    ncores = int(sys.argv[2])

datacards = [d for d in glob.glob(datacard_pattern) if "/cmb/" not in d]

cmds = ['card=DATACARD; echo $card; echo ${card/combined.txt.cmb/stage0.root}; echo ${card/combined.txt.cmb/prefit_shapes.root}; PostFitShapesFromWorkspace -w ${card/combined.txt.cmb/stage0.root}  -o ${card/combined.txt.cmb/prefit_shapes.root} -d ${card}'.replace("DATACARD",d) for d in datacards]

p = Pool(ncores)
p.map(execute, cmds)
