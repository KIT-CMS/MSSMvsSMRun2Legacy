#!/usr/bin/env python
# file is copied from old repository of HIG-17-020:
# https://raw.githubusercontent.com/KIT-CMS/CombineHarvester/SMHTTLegacy-dev/MSSMFull2016/scripts/plotMultiDimFit.py
import CombineHarvester.CombineTools.plotting as plot
import ROOT
import math
import argparse
from array import array

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs="+", help='Input files')
parser.add_argument(
    '--output', '-o', default='limit', help="""Name of the output
    plot without file extension""")
parser.add_argument(
    '--sm-exp', default=[], nargs="+", help="""Input files for the SM expectation""")
parser.add_argument(
    '--bg-exp', default=[], nargs="+", help="""Input files for the backgroud only expectation""")
parser.add_argument(
    '--cms-sub', default='Internal', help="""Text below the CMS logo""")
parser.add_argument(
    '--mass', default='', help="""Mass label on the plot""")
parser.add_argument(
    '--title-right', default='', help="""Right header text above the frame""")
parser.add_argument(
    '--title-left', default='', help="""Left header text above the frame""")
parser.add_argument(
    '--x-title', default='#sigma#font[42]{(gg#phi)}#upoint#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)', help="""Title for the x-axis""")
parser.add_argument(
    '--y-title', default='#sigma#font[42]{(bb#phi)}#upoint#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)', help="""Title for the x-axis""")
parser.add_argument(
    '--y-axis-max', default=None, help="""y-axis max""")
parser.add_argument(
    '--x-axis-max', default=None, help="""x-axis max""")
parser.add_argument(
    '--likelihood-database', action='store_true', help="""Output likelihood database instead of plot""")
parser.add_argument(
    '--debug-output', '-d', help="""If specified, write the contour TH2s and
    TGraphs into this output ROOT file""")
parser.add_argument(
    '--add-3sigma-contour', action='store_true', help="""Add 3 sigma contour to plots of likelihood
    scans.""")
parser.add_argument(
    "--interpolate-missing", action="store_true", help="""Interpolate missing values in the scan.""")
args = parser.parse_args()

#Create canvas and TH2D for each component
plot.ModTDRStyle(width=600, l=0.13, r=0.08)
ROOT.gStyle.SetNdivisions(510, 'XYZ')
plot.SetBirdPalette()
canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()

if args.debug_output is not None:
    debug = ROOT.TFile(args.debug_output, 'RECREATE')
else:
    debug = None

limit = plot.MakeTChain(args.files, 'limit')

graph = plot.TGraph2DFromTree(
    limit, "r_ggH", "r_bbH", '2*deltaNLL', 'quantileExpected > -0.5 && deltaNLL > 0 && deltaNLL < 1000')
best = plot.TGraphFromTree(
    limit, "r_ggH", "r_bbH", 'deltaNLL == 0')
plot.RemoveGraphXDuplicates(best)
# hists = plot.TH2FromTGraph2D(graph, method='BinCenterAligned')
hists = plot.TH2FromTGraph2D(graph, method='BinCenterAligned')
plot.fastFillTH2(hists, graph,interpolateMissing=args.interpolate_missing)
if args.bg_exp:
    limit_bg = plot.MakeTChain(args.bg_exp, 'limit')
    best_bg = plot.TGraphFromTree(
       limit_bg, "r_ggH", "r_bbH", 'deltaNLL == 0')
    plot.RemoveGraphXDuplicates(best_bg)
if args.sm_exp:
    limit_sm = plot.MakeTChain(args.sm_exp, 'limit')
    best_sm = plot.TGraphFromTree(
        limit_sm, "r_ggH", "r_bbH", 'deltaNLL == 0')
    plot.RemoveGraphXDuplicates(best_sm)
hists.SetMaximum(6)
hists.SetMinimum(0)
hists.SetContour(255)
# c2=ROOT.TCanvas()
# hists.Draw("COLZ")
# c2.SaveAs("heatmap.png")
#Set x and y axis maxima:
if args.y_axis_max is not None:
    y_axis_max = float(args.y_axis_max)
else:
    y_axis_max = float(hists.GetYaxis().GetXmax())

if args.x_axis_max is not None:
   x_axis_max = float(args.x_axis_max)
else:
   x_axis_max = float(hists.GetXaxis().GetXmax())
axis = ROOT.TH2D(hists.GetName(),hists.GetName(),hists.GetXaxis().GetNbins(),0,x_axis_max,hists.GetYaxis().GetNbins(),0,y_axis_max)
axis.Reset()
axis.GetXaxis().SetTitle(args.x_title)
axis.GetXaxis().SetLabelSize(0.035)
axis.GetYaxis().SetLabelSize(0.035)
axis.GetYaxis().SetTitle(args.y_title)
axis.GetXaxis().SetTitleSize(0.055)
axis.GetYaxis().SetTitleSize(0.055)
# Set ticks also on upper x and right y axis
for pad in pads:
    pad.SetTickx()
    pad.SetTicky()
ROOT.TGaxis.SetMaxDigits(3)
# ROOT.TGaxis.SetExponentOffset(-0.025, -0.05, "X")
# ROOT.TGaxis.SetExponentOffset(-0.037, -0.06, "X")

cont_1sigma = plot.contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.68, 2), 10, frameValue=20)
cont_2sigma = plot.contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.95, 2), 10, frameValue=20)
if args.add_3sigma_contour:
    cont_3sigma = plot.contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.997, 2), 10, frameValue=20)

if debug is not None:
    debug.WriteTObject(hists, 'hist')
    for i, cont in enumerate(cont_1sigma):
        debug.WriteTObject(cont, 'cont_1sigma_%i' % i)
    for i, cont in enumerate(cont_2sigma):
        debug.WriteTObject(cont, 'cont_2sigma_%i' % i)
    if args.add_3sigma_contour:
        for i, cont in enumerate(cont_3sigma):
            debug.WriteTObject(cont, 'cont_3sigma_%i' % i)

if args.sm_exp or args.bg_exp:
    legend = plot.PositionedLegend(0.5, 0.25, 3, 0.015)
else:
     legend = plot.PositionedLegend(0.3, 0.2, 3, 0.015)

legend.SetFillStyle(0)

pads[0].cd()
axis.Draw()
if args.add_3sigma_contour:
    for i, p in enumerate(cont_3sigma):
          p.SetLineStyle(1)
          p.SetLineWidth(2)
          p.SetLineColor(ROOT.kBlack)
          p.SetFillColor(ROOT.kCyan-10)
          p.SetFillStyle(1001)
          p.Draw("F SAME")
          p.Draw("L SAME")
          try:
              legend.AddEntry(cont_1sigma[0], "68% CL", "F")
          except IndexError:
              print("Problems with legend..")
              pass

for i, p in enumerate(cont_2sigma):
      p.SetLineStyle(1)
      p.SetLineWidth(2)
      p.SetLineColor(ROOT.kBlack)
      p.SetFillColor(ROOT.kBlue-10)
      p.SetFillStyle(1001)
      p.Draw("F SAME")
      p.Draw("L SAME")
      try:
          if args.add_3sigma_contour:
              legend.AddEntry(cont_2sigma[0], "95% CL", "F")
          else:
              legend.AddEntry(cont_1sigma[0], "68% CL", "F")
      except IndexError:
          print("Problems with legend..")
          pass

for i, p in enumerate(cont_1sigma):
      p.SetLineStyle(1)
      p.SetLineWidth(2)
      p.SetLineColor(ROOT.kBlack)
      p.SetFillColor(ROOT.kBlue-8)
      p.SetFillStyle(1001)
      p.Draw("F SAME")
      p.Draw("L SAME")
      try:
          if args.add_3sigma_contour:
              legend.AddEntry(cont_3sigma[0], "99.7% CL", "F")
          else:
              legend.AddEntry(cont_2sigma[0], "95% CL", "F")
      except IndexError:
          print("Problems with legend..")
          pass

best.SetMarkerStyle(34)
best.SetMarkerSize(3)
best.Draw("P SAME")
legend.AddEntry(best, "Best fit", "P")
if args.sm_exp:
    best_sm.SetMarkerStyle(33)
    best_sm.SetMarkerColor(1)
    best_sm.SetMarkerSize(3.0)
    best_sm.Draw("P SAME")
    legend.AddEntry(best_sm, "Expected for 125 GeV SM Higgs", "P")
if args.bg_exp:
    best_bg.SetMarkerStyle(33)
    best_bg.SetMarkerColor(46)
    best_bg.SetMarkerSize(3)
    best_bg.Draw("P SAME")
    legend.AddEntry(best_bg, "Expected for background only", "P")


if args.mass:
    legend.SetHeader("m_{#phi} = "+args.mass+" GeV")
legend.Draw("SAME")
if args.sm_exp:
    overlayLegend,overlayGraphs = plot.getOverlayMarkerAndLegend(legend, {legend.GetNRows()-1 : best_sm}, {legend.GetNRows()-1 : {"MarkerColor" : 2}}, markerStyle="P")

plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 11, 0.045, 0.035, 1.2, '', 1.0)
plot.DrawTitle(pads[0], args.title_right, 3)
plot.DrawTitle(pads[0], args.title_left, 1)
plot.FixOverlay()
if args.sm_exp:
    best_sm.Draw("P SAME")
    for overlayGraph in overlayGraphs:
        overlayGraph.Draw("P SAME")
    overlayLegend.Draw("SAME")
canv.Print('.pdf')
canv.Print('.png')
canv.Close()

if debug is not None:
    debug.Close()

if args.likelihood_database:
    output_file = open(args.output+".out","w")
    for i in range(0,limit.GetEntries()):
        limit.GetEntry(i)
        ggH = limit.r_ggH
        bbH = limit.r_bbH
        deltanll = limit.deltaNLL
        output_file.write("%(ggH)f %(bbH)f %(deltanll)f \n"%vars())
    output_file.close()
