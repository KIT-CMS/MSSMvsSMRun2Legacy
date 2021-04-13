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
        print "[WARNING] Command failed:", cmd


parser = argparse.ArgumentParser(
    description=
    "Script to run 'PostFitShapesFromWorkspace' in parallel for each category."
)
parser.add_argument(
    '--basedir',
    required=True,
    help="Path pattern to the folder containing all datacards of the analysis")
parser.add_argument('--workspacename',
                    default="ws.root",
                    help="Name of the workspace")
# parser.add_argument('--output_name', default="prefit_shapes.root", help = "Name of the output files")
parser.add_argument(
    '--freeze-arguments',
    default="",
    help=
    "Arguments to be frozen to a certain value. Use with options, e.g. '--freeze r=1'"
)
parser.add_argument(
    '--fit-arguments',
    default="",
    help=
    "Arguments to be used to apply fit result(s). Use with needed options, e.g. '-f <path-pattern-to-fitDiagnostics.root>:<fit-to-be-used> --sampling --postfit'"
)
parser.add_argument('--parallel',
                    type=int,
                    default=5,
                    help="Cores provided for parallel processing")
parser.add_argument('--dry_run',
                    action='store_true',
                    help="Don't execute, only list commands")

args = parser.parse_args()

datacards = [
    card for card in os.listdir(os.path.join(args.basedir, "combined", "cmb"))
    if "txt" in card and not "combined" in card
]
print("Running prefit shapes for {} histograms".format(len(datacards)))

basedir = args.basedir
cmds = []
for datacardfile in datacards:
    datacard = datacardfile.strip(".txt")
    category = datacard.split("_")[1]
    channel = datacard.split("_")[2]
    era = datacard.split("_")[3]
    workspace = os.path.join(basedir, era, datacard, args.workspacename)
    outputfile = os.path.join(basedir, "prefitshapes",
                              "{}_prefit_shape.root".format(datacard))
    if not os.path.exists(os.path.join(basedir, "prefitshapes")):
        os.makedirs(os.path.join(basedir, "prefitshapes"))
    restore_binningfile = os.path.join(basedir, "restore_binning",
                                       "{}.txt".format(datacard))
    commandstring = "PostFitShapesFromWorkspace -w {workspace} -o {outputfile} -d {restore_binningfile} {freezeargs} {fitargs}".format(
        workspace=workspace,
        outputfile=outputfile,
        restore_binningfile=restore_binningfile,
        freezeargs=args.freeze_arguments,
        fitargs=args.fit_arguments)
    cmds.append(commandstring)

p = Pool(args.parallel)
if args.dry_run:
    for cmd in cmds:
        print cmd
else:
    p.map(execute, cmds)
