#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import argparse
import numpy as np
import json
import os

parser = argparse.ArgumentParser(description="Compare two sets of hig-21-001 datacard shapes")
parser.add_argument("--first-shapefile", required=True, help="Path to the ROOT file with reference set of datacard shapes to compare with.")
parser.add_argument("--second-shapefile", required=True, help="Path to the ROOT file with a new set of datacard shapes.")
parser.add_argument("--consistency-2D-vs-1D", default=1, choices=[0,1], type=int, help="Flag to perform consistency check between 2D signal and 1D split categories. Default: %(default)s")
parser.add_argument("--check-uncertainties", default=1, choices=[0,1], type=int, help="Flag to perform check also uncertainties. Default: %(default)s")
parser.add_argument("--process", default=None, help="Compare specific process. Default: %(default)s")
parser.add_argument("--verbose", action="store_true", help="Flag to make verbose outputs")
parser.add_argument("--process-mapping", default=None, help="Path to a .json file with a name mapping from first to second shape file. The name stored in the first file is taken as reference name used as key. Default: %(default)s")

args = parser.parse_args()
ROOT.gROOT.SetBatch()

process_mapping = {}

if args.process_mapping:
    if not os.path.isfile(args.process_mapping):
        print("ERROR: file for process mapping does not exist. Aborting.")
        exit(1)
    process_mapping = json.load(open(args.process_mapping, 'r'))

first = ROOT.TFile.Open(args.first_shapefile,"read")
second = ROOT.TFile.Open(args.second_shapefile,"read")

first_categories = set([k.GetName() for k in first.GetListOfKeys()])
second_categories = set([k.GetName() for k in second.GetListOfKeys()])

common = first_categories.intersection(second_categories)

only_first = first_categories.difference(second_categories)
only_second = second_categories.difference(first_categories)

if len(only_first) > 0:
    print("Categories available only in "+args.first_shapefile)
    print(only_first)
if len(only_second) > 0:
    print("Categories available only in "+args.second_shapefile)
    print(only_second)

print("Checking now common categories")
for cat in common:
    print "\tConsidering",cat
    first_dir = first.Get(cat)
    second_dir = second.Get(cat)
    first_hists = set([k.GetName() for k in first_dir.GetListOfKeys()])
    second_hists = set([k.GetName() for k in second_dir.GetListOfKeys()])
    second_hists_mapped = set()
    for h in second_hists:
        for fst, snd in process_mapping.items():
            if snd in h:
                h = h.replace(snd, fst)
                break
        second_hists_mapped.add(h)
    second_hists = second_hists_mapped
    common_hists = first_hists.intersection(second_hists)
    only_first_hists = first_hists.difference(second_hists)
    only_second_hists = second_hists.difference(first_hists)
    common_processes = set([k for k in common_hists if not k.endswith("Up") and not k.endswith("Down")])
    only_first_processes = set([k for k in only_first_hists if not k.endswith("Up") and not k.endswith("Down")])
    only_second_processes = set([k for k in only_second_hists if not k.endswith("Up") and not k.endswith("Down")])

    if args.process:
        common_hists = set([k for k in common_hists if args.process in k])
        only_first_hists = set([k for k in only_first_hists if args.process in k])
        only_second_hists = set([k for k in only_second_hists if args.process in k])
        common_processes = set([k for k in common_processes if args.process in k])
        only_first_processes = set([k for k in only_first_processes if args.process in k])
        only_second_processes = set([k for k in only_second_processes if args.process in k])

    if not args.check_uncertainties:
        common_hists = common_processes
        only_first_hists = only_first_processes
        only_second_hists = only_second_processes

    if len(only_first_processes) > 0:
        print("\tProcesses available only in "+args.first_shapefile)
        print "\t",only_first_processes

    only_first_hists_single = only_first_hists
    for proc in only_first_processes:
        only_first_hists_single = only_first_hists_single.intersection(set([h for h in only_first_hists_single if not proc in h]))
    if len(only_first_hists_single) > 0:
        print("\tIndividual histograms available only in "+args.first_shapefile)
        print "\t",only_first_hists_single

    if len(only_second_processes) > 0:
        print("\tProcesses available only in "+args.second_shapefile)
        print "\t",only_second_processes

    only_second_hists_single = only_second_hists
    for proc in only_second_processes:
        only_second_hists_single = only_second_hists_single.intersection(set([h for h in only_second_hists_single if not proc in h]))
    if len(only_second_hists_single) > 0:
        print("\tIndividual histograms available only in "+args.second_shapefile)
        print "\t",only_second_hists_single

    for hist in common_hists:
        hist2 = hist
        for fst, snd in process_mapping.items():
            if fst in hist2:
                hist2 = str(hist2.replace(fst,snd))
                break
        first_hist = first_dir.Get(hist)
        second_hist = second_dir.Get(hist2)
        diff_hist = first_hist.Clone("diff_"+hist)
        diff_hist.Add(second_hist, -1.0)
        diff_vals = np.array([abs(diff_hist.GetBinContent(i+1)) for i in range(diff_hist.GetNbinsX())])
        diff_value = np.sum(diff_vals)
        if diff_value >  0.0:
            print "\t\tDifference spotted:",hist,diff_value,"ratio to first:",diff_value/first_hist.Integral()
            if args.verbose:
                print "\t\tDifference in bins:",hist,diff_vals

if args.consistency_2D_vs_1D:
    print("")
    first_sm_signal_categories = set([c for c in first_categories if "xxh" in c])
    second_sm_signal_categories = set([c for c in second_categories if "xxh" in c])

    # if more than 1, then split categories available for 2D xxh. Usually 1 "xxh" category per channel (ROOT files provided channel-wise)
    if len(first_sm_signal_categories) > 1:
        print("Checking consistency of 1D split SM signal categories with 2D SM signal category for",args.first_shapefile)
        assert(len([c for c in first_sm_signal_categories if "bin" not in c]) == 1)
        sm_signal_2D = [first.Get(c) for c in first_sm_signal_categories if "bin" not in c][0]
        sm_signal_1D = [first.Get(c) for c in sorted(first_sm_signal_categories) if "bin" in c]
        # There should be 6 splits
        assert(len(sm_signal_1D) == 6)
        for k in sm_signal_2D.GetListOfKeys():
            h_2D = sm_signal_2D.Get(k.GetName())
            sm_2D_values = np.array([h_2D.GetBinContent(i+1) for i in range(h_2D.GetNbinsX())])
            sm_1D_split_values = np.array([])
            for cat_1D in sm_signal_1D:
                h_1D = cat_1D.Get(k.GetName())
                sm_1D_split_values = np.concatenate((sm_1D_split_values, np.array([h_1D.GetBinContent(i+1) for i in range(h_1D.GetNbinsX())])))

            diff_2D_vs_split_1D = np.abs(sm_2D_values - sm_1D_split_values)
            if np.sum(diff_2D_vs_split_1D) > 0.0:
                print "\tDifference spotted:",k.GetName(),diff_2D_vs_split_1D

    # if more than 1, then split categories available for 2D xxh. Usually 1 "xxh" category per channel (ROOT files provided channel-wise)
    if len(second_sm_signal_categories) > 1:
        print("Checking consistency of 1D split SM signal categories with 2D SM signal category for "+args.second_shapefile)
        assert(len([c for c in second_sm_signal_categories if "bin" not in c]) == 1)
        sm_signal_2D = [second.Get(c) for c in second_sm_signal_categories if "bin" not in c][0]
        sm_signal_1D = [second.Get(c) for c in sorted(second_sm_signal_categories) if "bin" in c]
        # There should be 6 splits
        assert(len(sm_signal_1D) == 6)
        for k in sm_signal_2D.GetListOfKeys():
            h_2D = sm_signal_2D.Get(k.GetName())
            sm_2D_values = np.array([h_2D.GetBinContent(i+1) for i in range(h_2D.GetNbinsX())])
            sm_1D_split_values = np.array([])
            for cat_1D in sm_signal_1D:
                h_1D = cat_1D.Get(k.GetName())
                sm_1D_split_values = np.concatenate((sm_1D_split_values, np.array([h_1D.GetBinContent(i+1) for i in range(h_1D.GetNbinsX())])))

            diff_2D_vs_split_1D = np.abs(sm_2D_values - sm_1D_split_values)
            if np.sum(diff_2D_vs_split_1D) > 0.0:
                print "\tSplitting difference spotted:",k.GetName(),np.sum(diff_2D_vs_split_1D)
                print "\t2D      :",sm_2D_values
                print "\t1D split:",sm_1D_split_values
