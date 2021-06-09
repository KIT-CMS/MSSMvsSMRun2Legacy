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
parser.add_argument('--output-folder', required = True, help = "Main folder, where the datacards should be created")
parser.add_argument('--analysis', required = True, help = "Analysis to be prepared with morphing")
parser.add_argument('--sub-analysis', required = True, help = "Sub-analysis to be prepared with morphing")
parser.add_argument('--categorization', required = True, help = "Categorization type to be prepared with morphing")
parser.add_argument('--sm-like-hists', required = True, help = "Templates type for SM like Higss boson to be prepared with morphing")
parser.add_argument('--category-list', required = True, help = "Category list, which will be used for parallelization of morphing")
parser.add_argument('--variable', required = True, help = "Variable to be used for the list of categories")
parser.add_argument('--eras', required = True, help = "Eras list, which will be used for parallelization of morphing")
parser.add_argument('--parallel', type=int, default=5, help = "Cores provided for parallel morphing")
parser.add_argument('--additional-arguments', type=str, default="--auto_rebin=1" , help = "Additional arguments to be passed to the Morphing executable")
parser.add_argument('--dry-run',action='store_true', help = "Don't execute, only list Morphing commands")
parser.add_argument('--sm-gg-fractions',
                    default = '${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root',
                    help = "sm-gg-fractions file to use")
parser.add_argument('--sm',action='store_true', help = "If set to true, sm categories are used")
args = parser.parse_args()

categories = []
with open(args.category_list, "r") as f:
    categories = [l.strip() for l in f.readlines()]

eras = args.eras.split(',')

commands = []

command_template = "MorphingMSSMvsSM --era={ERA} --category={CATEGORY} --output_folder={OUTPUT}" \
                   " --analysis={ANALYSIS} --sub-analysis={SUB_ANALYSIS} --categorization={CATEGORIZATION} --sm-like-hists={SM_LIKE_HISTS}" \
                   " --variable={VARIABLE} --sm_gg_fractions={SM_GG_FRACTIONS} {ADDITIONALARGS}"

for era in eras:
    for category in categories:
        command = command_template.format(ERA=era, CATEGORY=category, ANALYSIS=args.analysis,
                                          ADDITIONALARGS=args.additional_arguments, OUTPUT=args.output_folder,
                                          VARIABLE=args.variable, SM_GG_FRACTIONS=args.sm_gg_fractions,
                                          SUB_ANALYSIS=args.sub_analysis, CATEGORIZATION=args.categorization, SM_LIKE_HISTS=args.sm_like_hists)
        commands.append(command)

if args.sm:
    commands = ["{} --sm=true".format(command) for command in commands]

if args.dry_run:
    for command in commands:
        print command

else:
    p = Pool(args.parallel)
    p.map(execute, commands)
