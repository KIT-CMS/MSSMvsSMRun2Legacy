import CombineHarvester.CombineTools.plotting as plot 
import ROOT
import math
import argparse
from array import array
import sys

ROOT.gROOT.SetBatch(ROOT.kTRUE)

sys.setrecursionlimit(20000)

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
    '--cms-sub', default='Supplementary', help="""Text below the CMS logo""")
parser.add_argument(
    '--mass', default='', help="""Mass label on the plot""")
parser.add_argument(
    '--title-right', default='138 fb^{-1} (13 TeV)', help="""Right header text above the frame""")
parser.add_argument(
    '--title-left', default='', help="""Left header text above the frame""")
parser.add_argument(
    '--x-title', default='g_{t}^{X} #sqrt{#font[52]{B}#font[42]{(X#rightarrow#tau#tau)}}', help="""Title for the x-axis""")
parser.add_argument(
    '--y-title', default='g_{b}^{X} #sqrt{#font[52]{B}#font[42]{(X#rightarrow#tau#tau)}}', help="""Title for the y-axis""")
parser.add_argument(
    '--x-min', default=None, help="minimum xaxis point")
parser.add_argument(
    '--x-max', default=None, help="maximum xaxis point")
parser.add_argument(
    '--y-min', default=None, help="minimum yaxis point")
parser.add_argument(
    '--y-max', default=None, help="maximum yaxis point")
parser.add_argument(
    '--debug-output', '-d', help="""If specified, write the contour TH2s and
    TGraphs into this output ROOT file""")
parser.add_argument('--kappa', action='store_true', help='Produce 2D plot for kappa parameterisation', default=False)
parser.add_argument('--CPodd', action='store_true', help='Produce scans of yukawa couplings', default=False)
args = parser.parse_args()

if args.CPodd: X='A'
else: X='H'

#Create canvas and TH2D for each component
plot.ModTDRStyle(width=600, l=0.12)
ROOT.gStyle.SetNdivisions(510, 'XYZ')
plot.SetBirdPalette()
canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()

if args.debug_output is not None:
    debug = ROOT.TFile(args.debug_output, 'RECREATE')
else:
    debug = None

limit = plot.MakeTChain(args.files, 'limit')
if args.CPodd:
  graph = plot.TGraph2DFromTree(
      limit, "Yt_A", "Yb_A", '2*deltaNLL', 'quantileExpected > -0.5 && deltaNLL > 0 && deltaNLL < 1000')
  best = plot.TGraphFromTree(
      limit, "Yt_A", "Yb_A", 'deltaNLL == 0')
else:
  graph = plot.TGraph2DFromTree(
      limit, "Yt_H", "Yb_H", '2*deltaNLL', 'quantileExpected > -0.5 && deltaNLL > 0 && deltaNLL < 1000')
  best = plot.TGraphFromTree(
      limit, "Yt_H", "Yb_H", 'deltaNLL == 0')
plot.RemoveGraphXDuplicates(best)
hists = plot.TH2FromTGraph2D(graph, method='BinCenterAligned')
plot.fastFillTH2(hists, graph,interpolateMissing=True)
if args.bg_exp:
    limit_bg = plot.MakeTChain(args.bg_exp, 'limit')
    best_bg = plot.TGraphFromTree(
        limit_bg, "alpha", "muF", 'deltaNLL == 0')
    plot.RemoveGraphXDuplicates(best_bg)

# If included just plot SM exp at 1,1
if args.sm_exp:
    #limit_sm = plot.MakeTChain(args.sm_exp, 'limit')
    #best_sm = ROOT.TGraph( 1, array('d', [1,]), array('d', [1,]))
    best_sm = ROOT.TGraph
    plot.RemoveGraphXDuplicates(best_sm)
hists.SetMaximum(6)
hists.SetMinimum(0)
hists.SetContour(255)
# c2=ROOT.TCanvas()
# hists.Draw("COLZ")
# c2.SaveAs("heatmap.png")

extra=0.

if int(args.mass) in [80,100,160,180,200]: extra=0.5
if int(args.mass) in [95]: extra=1.0

axis = ROOT.TH2D(hists.GetName(),hists.GetName(),hists.GetXaxis().GetNbins(),0,hists.GetXaxis().GetXmax()-0.1,hists.GetYaxis().GetNbins(),hists.GetYaxis().GetXmin(),hists.GetYaxis().GetXmax()+extra)
axis.Reset()
axis.GetXaxis().SetTitleOffset(0.9)
axis.GetXaxis().SetTitle(args.x_title.replace('X',X))
axis.GetXaxis().SetLabelSize(0.025)
axis.GetYaxis().SetLabelSize(0.025)
axis.GetYaxis().SetTitleOffset(0.9)
axis.GetYaxis().SetTitle(args.y_title.replace('X',X))

cont_1sigma = plot.contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.68, 2), 10, frameValue=20)
cont_2sigma = plot.contourFromTH2(hists, ROOT.Math.chisquared_quantile_c(1 - 0.95, 2), 10, frameValue=20)
print cont_1sigma
for item in cont_1sigma :
    print item
print cont_2sigma
for item in cont_2sigma :
    print item

if debug is not None:
    debug.WriteTObject(hists, 'hist')
    for i, cont in enumerate(cont_1sigma):
        debug.WriteTObject(cont, 'cont_1sigma_%i' % i)
    for i, cont in enumerate(cont_2sigma):
        debug.WriteTObject(cont, 'cont_2sigma_%i' % i)

if args.sm_exp or args.bg_exp:
    legend = plot.PositionedLegend(0.47, 0.20, 3, 0.015)
else:
    legend = plot.PositionedLegend(0.3, 0.2, 3, 0.015)

pads[0].cd()
pads[0].SetTicks(1)
axis.Draw()
for i, p in enumerate(cont_2sigma):
    p.SetLineStyle(1)
    p.SetLineWidth(2)
    p.SetLineColor(ROOT.kBlack)
    p.SetFillColor(ROOT.kBlue-10)
    p.SetFillStyle(1001)
    p.Draw("F SAME")
    p.Draw("L SAME")
    if i>0: continue
    legend.AddEntry(cont_1sigma[0], "68% CL", "F")

for i, p in enumerate(cont_1sigma):
    p.SetLineStyle(1)
    p.SetLineWidth(2)
    p.SetLineColor(ROOT.kBlack)
    p.SetFillColor(ROOT.kBlue-8)
    p.SetFillStyle(1001)
    p.Draw("F SAME")
    p.Draw("L SAME")
    if i>0: continue
    legend.AddEntry(cont_2sigma[0], "95% CL", "F")

best.SetMarkerStyle(34)
best.SetMarkerSize(3)
best.Draw("P SAME")
legend.AddEntry(best, "Best fit", "P")

if args.sm_exp:
    best_sm.SetMarkerStyle(33)
    best_sm.SetMarkerColor(1)
    best_sm.SetMarkerSize(3.0)
    best_sm.Draw("P SAME")
    legend.AddEntry(best_sm, "SM Expected", "P")
if args.bg_exp:
    best_bg.SetMarkerStyle(33)
    best_bg.SetMarkerColor(46)
    best_bg.SetMarkerSize(3)
    best_bg.Draw("P SAME")
    legend.AddEntry(best_bg, "Expected for background only", "P")


if args.mass:
    legend.SetHeader("m_{%(X)s} = " % vars() +args.mass+" GeV")
legend.Draw("SAME")
if args.sm_exp:
    overlayLegend,overlayGraphs = plot.getOverlayMarkerAndLegend(legend, {legend.GetNRows()-1 : best_sm}, {legend.GetNRows()-1 : {"MarkerColor" : 2}}, markerStyle="P")

plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 11, 0.045, 0.035, 1.2, '', 1.0)
pads[0].cd()
lumi = ROOT.TLatex(.7,1.05,"X fb^{-1} (13 TeV)")
lumi.SetTextSize(0.03)
#lumi.DrawTextNDC(.7,.96,"35.9 / fb (13 TeV)" )

plot.DrawTitle(pads[0], args.title_right, 3)
plot.DrawTitle(pads[0], args.title_left, 1)
plot.FixOverlay()
if args.sm_exp:
    best_sm.Draw("P SAME")
    for overlayGraph in overlayGraphs:
        print "test"
        overlayGraph.Draw("P SAME")
    overlayLegend.Draw("SAME")
canv.Print('.pdf')
canv.Print('.png')
canv.Print('.C')
canv.Close()

if debug is not None:
    debug.Close()
