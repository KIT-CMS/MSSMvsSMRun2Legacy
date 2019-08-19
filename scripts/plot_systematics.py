#!/usr/bin/env python
# -*- coding: utf-8 -*-


import ROOT as r
import sys
import os

fname = sys.argv[1]
category = sys.argv[2]
process = sys.argv[3]
outputfolder = sys.argv[4]

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)


f = r.TFile.Open(fname, "read")
d = f.Get(category)

systup = [ p.GetName() for p in d.GetListOfKeys() if process in p.GetName() and "Up" in p.GetName()]
systdown  = [ p.GetName() for p in d.GetListOfKeys() if process in p.GetName() and "Down" in p.GetName()]

if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

if not os.path.exists(os.path.join(outputfolder,category)):
    os.makedirs(os.path.join(outputfolder,category))

if not os.path.exists(os.path.join(outputfolder,category,process)):
    os.makedirs(os.path.join(outputfolder,category,process))

if len(systup) != len(systdown):
    print "Found different amount of Up & Down shapes for process %s in categoy %s. Up: %d, Down: %d. Aborting"%(process, category, len(systup), len(systdown))
    exit(1)

print "Making plots for %d shape systmatics."%len(systup)

c = r.TCanvas()
c.cd()

for s_u, s_d in zip(systup, systdown):
    c.Clear()
    nominal = d.Get(process)
    nominalcopy = d.Get(process).Clone()
    nominalcopy.SetFillColorAlpha(r.kGray,0.6)
    nominalcopy2 = d.Get(process).Clone()
    nominalcopy2.SetLineWidth(1)
    nominalcopy2.SetLineStyle(9)
    nominalcopy2.SetLineColor(r.kBlack)
    up = d.Get(s_u)
    up.SetLineWidth(3)
    up.SetLineColor(r.kRed+2)
    down = d.Get(s_d)
    down.SetLineWidth(2)
    down.SetLineColor(r.kCyan+2)
    up.Divide(nominal)
    down.Divide(nominal)
    nominalcopy.Divide(nominal)
    nominalcopy2.Divide(nominal)
    nominalcopy.GetYaxis().SetRangeUser(0.8,1.2)
    nominalcopy.Draw("][e2")
    nominalcopy2.Draw("][hist same")
    up.Draw("][hist same")
    down.Draw("][hist same")
    c.SaveAs(os.path.join(outputfolder,category,process,s_u.replace(process+"_","").replace("Up","")+".pdf"))
    c.SaveAs(os.path.join(outputfolder,category,process,s_u.replace(process+"_","").replace("Up","")+".png"))

f.Close()
