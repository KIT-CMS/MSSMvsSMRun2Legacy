#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit");
ROOT.gROOT.SetBatch()
import os
import CombineHarvester.CombineTools.ch as ch
import argparse

p = argparse.ArgumentParser("Script to print out systematics for a complete workspace")
p.add_argument("--workspace", required=True, help="Path to the workspace file")

args = p.parse_args()

f = ROOT.TFile.Open(args.workspace, "read")
ws = f.Get("w")

cmb = ch.CombineHarvester()
cmb.SetFlag("workspaces-use-clone", True)
ch.ParseCombineWorkspace(cmb, ws, "ModelConfig", "data_obs", False)
cmb.PrintParams()
