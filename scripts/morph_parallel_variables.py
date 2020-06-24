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
parser.add_argument('--eras', required = True, help = "Eras list, which will be used for parallelization of morphing")
parser.add_argument('--parallel', type=int, default=5, help = "Cores provided for parallel morphing")
parser.add_argument('--dry_run',action='store_true', help = "Don't execute, only list Morphing commands")

args = parser.parse_args()

variables_per_cat = {
    "inclusive" : {
        "mt" : ['mt_1_puppi','nbtag'],
        "et" : ['mt_1_puppi','nbtag'],
        "tt" : ['nbtag'],
        "em" : ['pZetaPuppiMissVis','nbtag'],
    },
    "signal_region" : {
        "mt" : ['nbtag'],
        "et" : ['nbtag'],
        "tt" : ['nbtag'],
        "em" : ['nbtag'],
    },
    "nobtag_lowmsv" : {
        "mt" : ['DiTauDeltaR','mjj','pt_tt_puppi','jdeta','njets'],
        "et" : ['DiTauDeltaR','mjj','pt_tt_puppi','jdeta','njets'],
        "tt" : ['DiTauDeltaR','mjj','pt_tt_puppi','jdeta','njets'],
        "em" : ['DiTauDeltaR','mjj','pt_tt_puppi','jdeta','njets'],
    },
}

eras = args.eras.split(',')

commands = []

command_template = "MorphingCatVariables --era={ERA} --category={CATEGORY}  --output_folder={OUTPUT} --variable={VARIABLE}"

for era in eras:
    for category in variables_per_cat:
            for channel in variables_per_cat[category]:
                for variable in variables_per_cat[category][channel]:
                    catname = "_".join([channel,category])
                    command = command_template.format(ERA=era, CATEGORY=catname, OUTPUT=args.output_folder, VARIABLE=variable)
                    commands.append(command)
                    print command

if not args.dry_run:
    p = Pool(args.parallel)
    p.map(execute, commands)
