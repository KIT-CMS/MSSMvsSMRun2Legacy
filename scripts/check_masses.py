#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
import argparse
import glob
import os

parser = argparse.ArgumentParser(
    description="Checking the mass ranges for a given folder of model files"
)

parser.add_argument(
    "--inputfolder",
    type=str,
    required=True,
    help="path to folder with model files following the 'mh*.root' ")

args = parser.parse_args()
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

flist = glob.glob(os.path.join(args.inputfolder,"m*13.root"))
higgses = ['h', 'A', 'H', 'H1', 'H2', 'H3']
for f in flist:
    F = ROOT.TFile.Open(f)
    scenario = os.path.basename(f).replace(".root","")
    print "Considering",scenario,"scenario"
    for higgs in higgses:
        histname = "m_"+higgs
        mhiggs = F.Get(histname)
        if mhiggs:
            print "\tRange for",histname,": [",mhiggs.GetMinimum(),",",mhiggs.GetMaximum(),"]"
        
