#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT as R
R.gROOT.SetBatch()
import sys

fpath = sys.argv[1]
poi1 = sys.argv[2]
poi2 = sys.argv[3]

F = R.TFile.Open(fpath,"read")
result = F.Get("fit_s")
params = result.floatParsInit()
correlation =  round(result.correlation(params.find(poi1),params.find(poi2)),3)
print "Determined correlation between {POI1} and {POI2}: {CORR}".format(POI1=poi1, POI2=poi2, CORR=correlation)
