#!/usr/bin/env python
import ROOT
import math
from functools import partial
import CombineHarvester.CombineTools.plotting as plot
import json
import argparse
import os.path

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)

plot.ModTDRStyle(width=700, l=0.13)
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetMarkerSize(0.7)

NAMECOUNTER = 0

def read(scan, param, files, ycut, rezero=True):
    goodfiles = [f for f in files if plot.TFileIsGood(f)]
    limit = plot.MakeTChain(goodfiles, 'limit')
    graph = plot.TGraphFromTree(limit, param, '2*deltaNLL', 'quantileExpected > -1.5')
    graph.SetName(scan)
    graph.Sort()
    plot.RemoveGraphXDuplicates(graph)
    plot.RemoveGraphYAbove(graph, ycut)
    if rezero: plot.ReZeroTGraph(graph, rezero)
    # graph.Print()
    return graph


def Eval(obj, x, params):
    return obj.Eval(x[0])


def BuildScan(scan, param, files, color, yvals, ycut):
    graph = read(scan, param, files, ycut)
    bestfit = None
    for i in xrange(graph.GetN()):
        if graph.GetY()[i] == 0.:
            bestfit = graph.GetX()[i]
    graph.SetMarkerColor(color)
    spline = ROOT.TSpline3("spline3", graph)
    global NAMECOUNTER
    func = ROOT.TF1('splinefn'+str(NAMECOUNTER), partial(Eval, spline), graph.GetX()[0], graph.GetX()[graph.GetN() - 1], 1)
    NAMECOUNTER += 1
    func.SetLineColor(color)
    func.SetLineWidth(3)
    assert(bestfit is not None)
    crossings = {}
    cross_1sig = None
    cross_2sig = None
    other_1sig = []
    other_2sig = []
    val = None
    val_2sig = None
    for yval in yvals:
        crossings[yval] = plot.FindCrossingsWithSpline(graph, func, yval)
        for cr in crossings[yval]:
            cr["contains_bf"] = cr["lo"] <= bestfit and cr["hi"] >= bestfit
    for cr in crossings[yvals[0]]:
        if cr['contains_bf']:
            val = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
            cross_1sig = cr
        else:
            other_1sig.append(cr)
    if len(yvals) > 1:
        for cr in crossings[yvals[1]]:
            if cr['contains_bf']:
                val_2sig = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
                cross_2sig = cr
            else:
                other_2sig.append(cr)
    else:
        val_2sig = (0., 0., 0.)
        cross_2sig = cross_1sig
    return {
        "graph"     : graph,
        "spline"    : spline,
        "func"      : func,
        "crossings" : crossings,
        "val"       : val,
        "val_2sig": val_2sig,
        "cross_1sig" : cross_1sig,
        "cross_2sig" : cross_2sig,
        "other_1sig" : other_1sig,
        "other_2sig" : other_2sig
    }

parser = argparse.ArgumentParser()

parser.add_argument('--obs', required=True, help='Observed input file for the scan')
parser.add_argument('--exp', required=True, help='Expected input file for the scan')
parser.add_argument('--y-cut', type=float, default=7., help='Remove points with y > y-cut')
parser.add_argument('--y-max', type=float, default=8., help='y-axis maximum')
parser.add_argument('--y-min', type=float, default=0., help='y-axis minimum')
parser.add_argument('--output', '-o', help='output name without file extension', default='scan')
parser.add_argument('--POI', help='use this parameter of interest', default='r')
parser.add_argument('--translate', default=None, help='json file with POI name translation')
parser.add_argument('--obs-label', default='Observed', type=str, help='legend label for the observed scan')
parser.add_argument('--obs-color', default=1, type=int, help='line and marker color for the observed scan')
parser.add_argument('--exp-label', default='Observed', type=str, help='legend label for the observed scan')
parser.add_argument('--exp-color', default=4, type=int, help='line and marker color for the observed scan')
parser.add_argument('--logo', default='CMS')
parser.add_argument('--logo-sub', default='Internal')
args = parser.parse_args()

print '--------------------------------------'
print  args.output
print '--------------------------------------'

fixed_name = args.POI
if args.translate is not None:
    with open(args.translate) as jsonfile:
        name_translate = json.load(jsonfile)
    if args.POI in name_translate:
        fixed_name = name_translate[args.POI]

yvals = [1., 4.]


obs_scan = BuildScan(args.output, args.POI, [args.obs], args.obs_color, yvals, args.y_cut)
exp_scan = BuildScan(args.output, args.POI, [args.exp], args.exp_color, yvals, args.y_cut)


canv = ROOT.TCanvas(args.output, args.output)
pads = plot.OnePad()

exp_scan['graph'].SetMarkerColor(args.exp_color)
exp_scan['graph'].SetMarkerStyle(ROOT.kCircle)
exp_scan['graph'].SetMarkerSize(1.)
exp_scan['graph'].SetLineWidth(3)
exp_scan['graph'].SetLineColor(args.exp_color)
#exp_scan['graph'].Draw('AP')
exp_scan['graph'].Draw('AL')

obs_scan['graph'].SetMarkerColor(args.obs_color)
obs_scan['graph'].SetMarkerStyle(ROOT.kCircle)
obs_scan['graph'].SetMarkerSize(1.)
obs_scan['graph'].SetLineWidth(3)
obs_scan['graph'].SetLineColor(args.obs_color)

axishist = plot.GetAxisHist(pads[0])

axishist.SetMinimum(args.y_min)
axishist.SetMaximum(args.y_max)
axishist.GetYaxis().SetTitle("#minus 2 #upoint #Delta ln L")
axishist.GetXaxis().SetTitle("%s" % fixed_name)
axishist.GetXaxis().SetTitleSize(0.07)
axishist.GetXaxis().SetTitleOffset(0.7)

initial_min = min([g.GetX()[0] for g in [obs_scan['graph'],exp_scan['graph']]])
initial_max = max([g.GetX()[g.GetN()-1] for g in [obs_scan['graph'],exp_scan['graph']]])
distance = initial_max - initial_min
new_min = initial_min - 0.05*distance
new_max = initial_max + 0.05*distance
axishist.GetXaxis().SetLimits(new_min,new_max)
axishist.Draw("axis")

line = ROOT.TLine()
line.SetLineStyle(2)
line.SetLineWidth(2)

for yval in yvals:
    line.SetLineColor(ROOT.kRed)
    plot.DrawHorizontalLine(pads[0], line, yval)
    line.SetLineColor(args.exp_color)
    for cr in exp_scan['crossings'][yval]:
        if cr['valid_lo']: line.DrawLine(cr['lo'], 0, cr['lo'], yval)
        if cr['valid_hi']: line.DrawLine(cr['hi'], 0, cr['hi'], yval)
    if args.obs:
        line.SetLineColor(ROOT.kRed)
        for cr in obs_scan['crossings'][yval]:
            if cr['valid_lo']: line.DrawLine(cr['lo'], 0, cr['lo'], yval)
            if cr['valid_hi']: line.DrawLine(cr['hi'], 0, cr['hi'], yval)


#exp_scan['graph'].Draw('LP same')
exp_scan['func'].Draw('l same')
obs_scan['graph'].Draw('LP same')
#obs_scan['graph'].Draw('L same')
#obs_scan['func'].Draw('l same')


box = ROOT.TBox(axishist.GetXaxis().GetXmin(), 0.625*args.y_max, axishist.GetXaxis().GetXmax(), args.y_max)
box.Draw()
pads[0].SetTopMargin(0.07)
pads[0].GetFrame().Draw()
pads[0].RedrawAxis()

crossings = obs_scan['crossings']
val_nom = obs_scan['val']
val_2sig = obs_scan['val_2sig']

exp_crossings = exp_scan['crossings']
exp_val_nom = exp_scan['val']
exp_val_2sig = exp_scan['val_2sig']

textfit = '%s = %.4f{}^{#plus %.4f}_{#minus %.4f}' % (fixed_name, val_nom[0], val_nom[1], abs(val_nom[2]))
exptextfit = '{}^{(#plus %.4f)}_{(#minus %.4f)}' % (exp_val_nom[1], abs(exp_val_nom[2]))


pt = ROOT.TPaveText(0.59, 0.82, 0.95, 0.88, 'NDCNB')
pt.AddText(textfit)
exppt = ROOT.TPaveText(0.78, 0.58, 0.95, 0.88, 'NDCNB')
exppt.AddText(exptextfit)

#pt.SetTextAlign(11)
pt.SetFillStyle(0)
pt.SetLineWidth(0)
pt.SetTextFont(42)
pt.SetTextSize(0.07)
pt.Draw()

exppt.SetFillStyle(0)
exppt.SetLineWidth(0)
exppt.SetTextFont(42)
exppt.SetTextColor(args.exp_color)
exppt.SetTextSize(0.07)
exppt.Draw()

plot.DrawCMSLogo(pads[0], args.logo, args.logo_sub, 11, 0.045, 0.035, 1.2,  cmsTextSize = 1.)
plot.DrawTitle(pads[0], "137 fb^{-1} (13 TeV)", 3)

legend = ROOT.TLegend(0.15, 0.62, 0.45, 0.76, '', 'NBNDC')
legend.SetFillStyle(0)

legend.AddEntry(obs_scan['func'], args.obs_label, 'L')
legend.AddEntry(exp_scan['func'], args.exp_label, 'L')
legend.Draw()

save_graph = obs_scan['graph'].Clone()
save_graph.GetXaxis().SetTitle('%s = %.3f %+.3f/%+.3f' % (fixed_name, val_nom[0], val_nom[2], val_nom[1]))
exp_graph = exp_scan['graph'].Clone()
exp_graph.GetXaxis().SetTitle('%s = %.3f %+.3f/%+.3f' % (fixed_name, exp_val_nom[0], exp_val_nom[2], exp_val_nom[1]))
outfile = ROOT.TFile(args.output+'.root', 'RECREATE')
outfile.WriteTObject(save_graph)
outfile.WriteTObject(exp_graph)
outfile.Close()
canv.Print('.pdf')
canv.Print('.png')

