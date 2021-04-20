#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob
import sys
import argparse
from multiprocessing import Pool

def execute(cmd):
    try:
        os.system(cmd)
    except:
        print "[WARNING] Command failed:",cmd

parser = argparse.ArgumentParser( description = "Script to run 'PostFitShapesFromWorkspace' in parallel for each category.")
parser.add_argument('--datacard_pattern', required = True, help = "Path pattern to the 'combined.txt.cmb' datacards to be used for the histograms")
parser.add_argument('--workspace_name', required = True, help = "Workspace name to be used to derive systematic uncertainties")
parser.add_argument('--output_name', default="prefit_shapes.root", help = "Name of the output files")
parser.add_argument('--freeze_arguments', default="", help = "Arguments to be frozen to a certain value. Use with options, e.g. '--freeze r=1'")
parser.add_argument('--fit_arguments', default="",
                    help = "Arguments to be used to apply fit result(s). Use with needed options, e.g. '-f <path-pattern-to-fitDiagnostics.root>:<fit-to-be-used> --sampling --postfit'")
parser.add_argument('--parallel', type=int, default=5, help = "Cores provided for parallel processing")
parser.add_argument('--dry_run',action='store_true', help = "Don't execute, only list commands")

args = parser.parse_args()

datacards = [d for d in glob.glob(args.datacard_pattern) if "/cmb/" not in d]

cmds = ['card=DATACARD; basedir=$(dirname $(dirname ${card})); maindir=$(dirname ${basedir}) category=$(basename $(dirname ${card})); echo $card; echo ${card/combined.txt.cmb/WSNAME}; echo ${card/combined.txt.cmb/OUTPUT}; PostFitShapesFromWorkspace -w ${card/combined.txt.cmb/WSNAME}  -o ${card/combined.txt.cmb/OUTPUT} -d ${maindir}/restore_binning/${category}.txt FREEZEARGS FITARGS'.replace("DATACARD",d).replace("FREEZEARGS",args.freeze_arguments).replace("FITARGS",args.fit_arguments).replace("WSNAME",args.workspace_name).replace("OUTPUT",args.output_name) for d in datacards]

p = Pool(args.parallel)
if args.dry_run:
    for cmd in cmds:
        print cmd
else:
    p.map(execute, cmds)
