#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import sys
import glob
import argparse
from multiprocessing import Pool

def execute(cmd):
    try:
        os.system(cmd)
    except:
        print "[WARNING] Command failed:",cmd

parser = argparse.ArgumentParser( description = "Compare Integrals of Processes between ML and Cutbased shapes")
parser.add_argument('--output_folder', required = True, help = "Main folder, where the datacards should be created")
parser.add_argument('--analysis', required = True, help = "Analysis to be prepared with morphing")
parser.add_argument('--category_list', required = True, help = "Category list, which will be used for parallelization of morphing")
parser.add_argument('--variable', required = True, help = "Variable to be used for the list of categories")
parser.add_argument('--eras', required = True, help = "Eras list, which will be used for parallelization of morphing")
parser.add_argument('--parallel', type=int, default=5, help = "Cores provided for parallel morphing")
parser.add_argument('--additional_arguments', type=str, default="--auto_rebin=1" , help = "Additional arguments to be passed to the Morphing executable")
parser.add_argument('--dry_run',action='store_true', help = "Don't execute, only list Morphing commands")
parser.add_argument('--sm_gg_fractions',
                    default = '${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3_mssm_mode.root',
                    help = "sm_gg_fractions file to use")
parser.add_argument('--sm',action='store_true', help = "If set to true, sm categories are used")
args = parser.parse_args()

categories = []
with open(args.category_list, "r") as f:
    categories = [l.strip() for l in f.readlines()]

eras = args.eras.split(',')

commands = []

command_template = "MorphingMSSMvsSM --era={ERA} --category={CATEGORY} --analysis={ANALYSIS} {ADDITIONALARGS} --output_folder={OUTPUT} --variable={VARIABLE} --sm_gg_fractions={SM_GG_FRACTIONS}"

for era in eras:
    for category in categories:
        command = command_template.format(ERA=era, CATEGORY=category, ANALYSIS=args.analysis, ADDITIONALARGS=args.additional_arguments, OUTPUT=args.output_folder, VARIABLE=args.variable, SM_GG_FRACTIONS=args.sm_gg_fractions)
        commands.append(command)

if args.sm:
    commands = ["{} --sm=true".format(command) for command in commands]

if args.dry_run:
    for command in commands:
        print command

else:
    p = Pool(args.parallel)
    p.map(execute, commands)
