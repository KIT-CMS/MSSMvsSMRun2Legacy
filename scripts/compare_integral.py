#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import ROOT as r
import sys
import os
import argparse
from prettytable import PrettyTable

parser = argparse.ArgumentParser(
    description="Compare Integrals of Processes between ML and Cutbased shapes"
)

parser.add_argument(
    "--file1",
    type=str,
    default=
    "/portal/ekpbms1/home/jbechtel/postprocessing/sync/sm-htt-analysis/htt_et.inputs-sm-Run2017-ML.root",
    help="path to first shapes")
parser.add_argument(
    "--file2",
    type=str,
    default=
    "/portal/ekpbms1/home/sbrommer/shape-producer/2016/sm-htt-analysis/CMSSW_10_2_13/src/CombineHarvester/MSSMvsSMRun2Legacy/shapes/htt_et.inputs-mssm-vs-sm-Run2017-m_sv_puppi.root",
    help="path to second shapes")
parser.add_argument("--type1",
                    type=str,
                    default="ml",
                    help="type of first shapes")
parser.add_argument("--type2",
                    type=str,
                    default="cut",
                    help="type of second shapes")
parser.add_argument('--channel',
                    type=str,
                    default="mt",
                    help="Channel to analyze")
parser.add_argument('--nominal',
                    action='store_true',
                    help=" Set to only compare nominals")
parser.add_argument('--diff',
                    action='store_true',
                    help=" Only print differences")
args = parser.parse_args()
r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

category_dict = {
    'ml': {
        'mt': ["ggh", "qqh", "ztt", "zll", "w", "tt", "ss", "misc"],
        'et': ["ggh", "qqh", "ztt", "zll", "w", "tt", "ss", "misc"],
        'tt': ["ggh", "qqh", "ztt", "noniso", "misc"],
        'em': ["ggh", "qqh", "ztt", "tt", "ss", "misc", "db"]
    },
    'cut': {
        'mt': ['nobtag', 'btag', 'wjets_control'],
        'et': ['nobtag', 'btag', 'wjets_control'],
        'tt': ['nobtag', 'btag'],
        'em': ['nobtag', 'btag', 'ttbar_control']
    }
}


def get_integrals(channel, file, nominal, type):
    categories = [
        '{}_{}'.format(channel, cat) for cat in category_dict[type][channel]
    ]
    rfile = r.TFile.Open(file, "read")
    if nominal:
        processes = [
            p.GetName() for p in rfile.Get(categories[0]).GetListOfKeys()
            if "Up" not in p.GetName() and "Down" not in p.GetName()
        ]
    else:
        processes = [
            p.GetName() for p in rfile.Get(categories[0]).GetListOfKeys()
        ]
    results = {}
    for i, process in enumerate(processes):
        print(' Reading {} Process {} /{}'.format(type, i + 1, len(processes)),
              end='\r')
        integral = 0
        for category in categories:
            d = rfile.Get(category)
            integral += d.Get(process).Integral()
        results[process] = integral
    rfile.Close()
    return results


print("Comparing Integral of two synced Shapes")
print("File1 : {} - Type {}".format(args.file1,args.type1))
print("File2 : {} - Type {}".format(args.file2,args.type2))

channel = args.channel
nominal = args.nominal
diff = args.diff
file1 = get_integrals(channel, args.file1, nominal, args.type1)
file2 = get_integrals(channel, args.file2, nominal, args.type2)

print(' --- results for {}---'.format(channel))
t = PrettyTable([
    'Process', 'File 1 Integral ({})'.format(args.type1),
    ' File 2 Integral ({})'.format(args.type2), 'Ratio'
])
for category in list(set(file1.keys()).intersection(file2.keys())):
    if file1[category] == 0 or file2[category] == 0:
        t.add_row([
            category,
            round(file1[category], 3),
            round(file2[category], 3), 0
        ])
    else:
        if diff:
            if round(file1[category] / file2[category], 3) != 1:
                t.add_row([
                    category,
                    round(file1[category], 3),
                    round(file2[category], 3),
                    round(file1[category] / file2[category], 3)
                ])
        else:
            t.add_row([
                category,
                round(file1[category], 3),
                round(file2[category], 3),
                round(file1[category] / file2[category], 3)
            ])
print(t.get_string(sortby='Ratio'))