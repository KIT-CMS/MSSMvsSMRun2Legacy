#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit");
ROOT.gROOT.SetBatch()
import os
import CombineHarvester.CombineTools.ch as ch
import argparse

p = argparse.ArgumentParser("Script to print out systematics of a given process in a given category for a given 'combined.txt.cmb' datacard. This does not included additional uncertainties, that could be added by the signal model used (e.g. by MSSM models)")
p.add_argument("--category", required=True, help="Category to be inspected")
p.add_argument("--process", required=True, help="Process to be inspected")
p.add_argument("--datacard", required=True, help="Path to the datacard 'combined.txt.cmb'")

args = p.parse_args()

cmb_card = ch.CombineHarvester()
cmb_card.SetFlag("workspaces-use-clone", True)
cmb_card.ParseDatacard(args.datacard, "", "", "", 0, "200")
cmb_card.bin([args.category]).cp().PrintProcs()
cmb_card.bin([args.category]).cp().process([args.process]).PrintSysts()
