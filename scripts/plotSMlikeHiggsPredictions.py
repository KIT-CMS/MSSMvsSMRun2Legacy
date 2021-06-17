#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT as R
#import CombineHarvester.CombineTools.plotting as plot
import plotting as plot
import numpy as np
import json
import argparse
import os
from array import array
R.gROOT.SetBatch()

plot.ModTDRStyle()
R.TColor.CreateGradientColorTable(3,
                                  array("d",[0.  ,0.5,  1.]),
                                  array("d",[0.  , 1.,  1.]),
                                  array("d",[0.35, 1.,0.65]),
                                  array("d",[1.  , 1.,  0.]),
                                  1000,  1.0)

parser = argparse.ArgumentParser(description="Derive comparisons for predictions for SM-like Higgs boson")
parser.add_argument('--mssm-benchmark', required=True, help="Path to the MSSM ROOT file for the benchmark scenario to be tested")
parser.add_argument('--bsm-sm-like', required=True, choices=["h", "H", "H1"], help="Name of the Higgs boson in the MSSM, which is supposed to be SM-like.")
parser.add_argument('--sm-predictions', required=True, help="Path to the .json file containing SM predictions for the SM-like Higgs boson.")
parser.add_argument('--plots', default="plots", help="Output directory for plots. Default: %(default)s")

args = parser.parse_args()

bsm_model = R.TFile.Open(args.mssm_benchmark, "read")
bsm_name = os.path.basename(args.mssm_benchmark.strip(".root"))

print("MODEL:",bsm_name)

sm_predictions = {}
with open(args.sm_predictions, "r") as smf:
    sm_predictions = json.load(smf)

C = R.TCanvas()
C.SetLeftMargin(1.3)
C.SetRightMargin(0.2)
C.cd()

bsm_predictions = {}
contour_graphs = {}

shift = 1e-7
mass_borders = [122.0,128.0]

bsm_model_names = {
    "mh125_13" : "M_{h}^{125}",
    "mh125EFT_13" : "M_{h,EFT}^{125}",
    "mh125_lc_13" : "M_{h}^{125}(#tilde{#chi})",
    "mh125EFT_lc_13" : "M_{h,EFT}^{125}(#tilde{#chi})",
    "mh125_ls_13" : "M_{h}^{125}(#tilde{#tau})",
    "mh125_align_13" : "M_{h}^{125}(alignment)",
    "mHH125_13" : "M_{H}^{125}(alignment)",
    "mh1125_CPV_13" : "M_{h_{1}}^{125}(CPV)",
    "mh125_muneg_1_13" : "M_{h}^{125}(#mu = #minus1 TeV)",
    "mh125_muneg_2_13" : "M_{h}^{125}(#mu = #minus2 TeV)",
    "mh125_muneg_3_13" : "M_{h}^{125}(#mu = #minus3 TeV)",
}

sf_range = [0.9, 1.1]
sf_contours = {0.9 : R.kBlue, 0.95 : R.kViolet-6, 0.99 : R.kCyan+1, 1.0 : R.kGreen+2, 1.02 : R.kRed+1, 1.1 : R.kMagenta, 1.3 : R.kBlack}
r_sf_contours = {0.98 : R.kBlue, 1.0 : R.kViolet-6, 1.02 : R.kCyan+1, 1.05 : R.kGreen+2, 1.1 : R.kRed+1, 1.3 : R.kMagenta, 1.5 : R.kBlack}

quantity_settings = {
    "sf_gg_{PHI}" : {
        "range" : sf_range,
        "contours" : sf_contours,
        "name" : "SF(gg#rightarrow{PHI}#rightarrow#tau#tau)".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "sf_bb_{PHI}" : {
        "range" : sf_range,
        "contours" : sf_contours,
        "name" : "SF(bb#rightarrow{PHI}#rightarrow#tau#tau)".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "sf_qq_{PHI}" : {
        "range" : sf_range,
        "contours" : sf_contours,
        "name" : "SF(qq#rightarrow{PHI}#rightarrow#tau#tau)".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "sf_gg_{PHI}_mass-corr" : {
        "range" : sf_range,
        "contours" : sf_contours,
        "name" : "mass-corrected SF(gg#rightarrow{PHI}#rightarrow#tau#tau)".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "sf_qq_{PHI}_mass-corr" : {
        "range" : sf_range,
        "contours" : sf_contours,
        "name" : "mass-corrected SF(qq#rightarrow{PHI}#rightarrow#tau#tau)".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "sf_bb_{PHI}_mass-corr" : {
        "range" : sf_range,
        "contours" : sf_contours,
        "name" : "mass-corrected SF(bb#rightarrow{PHI}#rightarrow#tau#tau)".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "r_sf_gg_{PHI}" : {
        "range" : sf_range,
        "contours" : r_sf_contours,
        "name" : "SF ratio for gg#rightarrow{PHI}#rightarrow#tau#tau".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "r_sf_qq_{PHI}" : {
        "range" : sf_range,
        "contours" : r_sf_contours,
        "name" : "SF ratio for qq#rightarrow{PHI}#rightarrow#tau#tau".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
    "r_sf_bb_{PHI}" : {
        "range" : sf_range,
        "contours" : r_sf_contours,
        "name" : "SF ratio for bb#rightarrow{PHI}#rightarrow#tau#tau".format(PHI="h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
    },
}

quantities = ['br_{PHI}_tautau', 'width_{PHI}', 'xs_gg_{PHI}', 'xs_bb_{PHI}']
quantities_for_plotting = []

quantity_range = [0.7, 1.3]
quantity_contours = {0.7 : R.kBlue, 0.8 : R.kViolet-6, 0.9 : R.kCyan+1, 1.0 : R.kGreen+2, 1.1 : R.kRed+1, 1.2 : R.kMagenta, 1.3 : R.kBlack}

quantity_settings['br_{PHI}_tautau_non-mass'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "BSM contributions to BR({PHI}#rightarrow#tau#tau)".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['br_{PHI}_tautau_mass-only'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "Mass effects on BR({PHI}#rightarrow#tau#tau)".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['xs_gg_{PHI}_non-mass'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "BSM contributions to #sigma(gg#rightarrow{PHI})".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['xs_gg_{PHI}_mass-only'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "Mass effects on #sigma(gg#rightarrow{PHI})".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['xs_bb_{PHI}_non-mass'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "BSM contributions to #sigma(bb#rightarrow{PHI})".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['xs_bb_{PHI}_mass-only'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "Mass effects on #sigma(bb#rightarrow{PHI})".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['width_{PHI}_non-mass'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "BSM contributions to #Gamma_{{PHI}}^{tot}".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['width_{PHI}_mass-only'] = {
    "range" : quantity_range,
    "contours" : quantity_contours,
    "name" : "Mass effects on #Gamma_{{PHI}}^{tot}".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

quantity_settings['m_{PHI}'] = {
    "range" : mass_borders,
    "contours" : {110.0 : R.kBlue, 120.5  : R.kViolet-6, 122.0 : R.kCyan+1, 125.0 : R.kGreen+2, 128.0 : R.kRed+1, 129.5 : R.kMagenta, 140.0 : R.kBlack},
    "name" : "m_{{PHI}}".replace("{PHI}","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),
}

if args.bsm_sm_like in ['h']:
    quantity_settings['gsq_{PHI}_VV'] = {
        "range" : [0.99, 1.01],
        "contours" : {0.97 : R.kBlue, 0.98 : R.kViolet-6, 0.99 : R.kCyan+1, 1.0 : R.kGreen+2, 1.01 : R.kRed+1, 1.02 : R.kMagenta, 1.04 : R.kBlack},
        "name" : "sin^{2}(#beta #minus #alpha)",
    }

if args.bsm_sm_like in ['H']:
    quantity_settings['gsq_{PHI}_VV'] = {
        "range" : [0.9, 0.92],
        "contours" : {0.88 : R.kBlue, 0.89 : R.kViolet-6, 0.90 : R.kCyan+1, 0.91 : R.kGreen+2, 0.92 : R.kRed+1, 0.93 : R.kMagenta, 0.94 : R.kBlack},
        "name" : "cos^{2}(#beta #minus #alpha)",
    }

for quantity in quantities:
    for postfix in "", "_SM":
        bsm_predictions[quantity + postfix] = bsm_model.Get((quantity + postfix).format(PHI=args.bsm_sm_like))

# Currently, these entries not available in the input file
if "EFT" in bsm_name or "mHH" in bsm_name:
    bsm_predictions['br_{PHI}_tautau_SM'] = bsm_predictions['br_{PHI}_tautau'].Clone('br_{PHI}_tautau_SM'.format(PHI=args.bsm_sm_like))
    bsm_predictions['width_{PHI}_SM'] = bsm_predictions['width_{PHI}'].Clone('width_{PHI}_SM'.format(PHI=args.bsm_sm_like))

# Currently, this entry is buggy in the input file
if "h1125" in bsm_name:
    bsm_predictions['xs_gg_{PHI}_SM'] = bsm_predictions['xs_gg_{PHI}'].Clone('xs_gg_{PHI}_SM'.format(PHI=args.bsm_sm_like))

bsm_predictions["m_{PHI}"] = bsm_model.Get("m_{PHI}".format(PHI=args.bsm_sm_like))
bsm_predictions["m_{PHI}_inverted"] = bsm_predictions["m_{PHI}"].Clone("m_{PHI}_inverted".format(PHI=args.bsm_sm_like))

if args.bsm_sm_like in ['h', 'H']:
    bsm_predictions["rescale_gt_H"] = bsm_model.Get("rescale_gt_H")

for bsm_pred in bsm_predictions.values():
    bsm_pred.SetContour(1000)

NXBins = bsm_predictions["m_{PHI}"].GetNbinsX() # mA or mHp
NYBins = bsm_predictions["m_{PHI}"].GetNbinsY() # tanb

# Transform m_{PHI} inverted histogram
for i_X in range(1, NXBins+1):
    for i_Y in range(1,NYBins+1):
        value = 0
        if bsm_predictions["m_{PHI}_inverted"].GetBinContent(i_X,i_Y) != 0:
            value = 1.-(1./bsm_predictions["m_{PHI}_inverted"].GetBinContent(i_X,i_Y))
        bsm_predictions["m_{PHI}_inverted"].SetBinContent(i_X,i_Y,value)

# Computing weight for VV coupling of SM-like Higgs (for CP conserving scenarios)

if args.bsm_sm_like in ['h', 'H']:
    bsm_predictions["gsq_{PHI}_VV"] = bsm_predictions["rescale_gt_H"].Clone("gsq_{PHI}_VV".format(PHI=args.bsm_sm_like))
    for i_X in range(1, NXBins+1):
        for i_Y in range(1,NYBins+1):
            tanb = bsm_predictions["gsq_{PHI}_VV"].GetYaxis().GetBinLowEdge(i_Y)
            yt_H = bsm_predictions["rescale_gt_H"].GetBinContent(i_X, i_Y)
            beta = np.arctan(tanb)
            alpha = np.arcsin(min(1.,max(-1.,yt_H*np.sin(beta))))
            if args.bsm_sm_like == "h":
                bsm_predictions["gsq_{PHI}_VV"].SetBinContent(i_X,i_Y,np.sin(beta - alpha)**2)
            elif args.bsm_sm_like == "H":
                bsm_predictions["gsq_{PHI}_VV"].SetBinContent(i_X,i_Y,np.cos(beta - alpha)**2)

# Compute non-mass contributions to quantities by dividing their SM-like equivalents at correct mass

for quantity in quantities:
    keyname = quantity + "_non-mass"
    bsm_predictions[keyname] = bsm_predictions[quantity].Clone(keyname.format(PHI=args.bsm_sm_like))
    bsm_predictions[keyname].Divide(bsm_predictions[quantity + "_SM"])
    
# Compute mass-only contributions to quantities by dividing SM-like quantities by the ones for SMH125
for quantity in quantities:
    keyname = quantity + "_mass-only"
    bsm_predictions[keyname] = bsm_predictions[quantity + "_SM"].Clone(keyname.format(PHI=args.bsm_sm_like))
    bsm_predictions[keyname].Scale(1. / sm_predictions[quantity.format(PHI="SMH125")])

# Compute total scale factors for ggPHI and qqPHI without mass correction (assuming signal templates are scaled to SMH125)
sfname = "sf_gg_{PHI}"
bsm_predictions[sfname] = bsm_predictions["br_{PHI}_tautau"].Clone(sfname.format(PHI=args.bsm_sm_like))
bsm_predictions[sfname].Multiply(bsm_predictions["xs_gg_{PHI}"])
bsm_predictions[sfname].Scale(1. / (sm_predictions["br_{PHI}_tautau".format(PHI="SMH125")] * sm_predictions["xs_gg_{PHI}".format(PHI="SMH125")]) )

sfname = "sf_bb_{PHI}"
bsm_predictions[sfname] = bsm_predictions["br_{PHI}_tautau"].Clone(sfname.format(PHI=args.bsm_sm_like))
bsm_predictions[sfname].Multiply(bsm_predictions["xs_bb_{PHI}"])
bsm_predictions[sfname].Scale(1. / (sm_predictions["br_{PHI}_tautau".format(PHI="SMH125")] * sm_predictions["xs_bb_{PHI}".format(PHI="SMH125")]) )

sfname = "sf_qq_{PHI}"
bsm_predictions[sfname] = bsm_predictions["br_{PHI}_tautau"].Clone(sfname.format(PHI=args.bsm_sm_like))
if args.bsm_sm_like in ['h', 'H']:
    bsm_predictions[sfname].Multiply(bsm_predictions["gsq_{PHI}_VV"])
bsm_predictions[sfname].Scale(1. / sm_predictions["br_{PHI}_tautau".format(PHI="SMH125")])

# Compute total scale factors for ggPHI and qqPHI with mass correction (assuming signal templates are scaled to SMH125)
sfname = "sf_gg_{PHI}_mass-corr"
bsm_predictions[sfname] = bsm_predictions["br_{PHI}_tautau_non-mass"].Clone(sfname.format(PHI=args.bsm_sm_like))
bsm_predictions[sfname].Multiply(bsm_predictions["xs_gg_{PHI}_non-mass"])

sfname = "sf_bb_{PHI}_mass-corr"
bsm_predictions[sfname] = bsm_predictions["br_{PHI}_tautau_non-mass"].Clone(sfname.format(PHI=args.bsm_sm_like))
bsm_predictions[sfname].Multiply(bsm_predictions["xs_bb_{PHI}_non-mass"])

sfname = "sf_qq_{PHI}_mass-corr"
bsm_predictions[sfname] = bsm_predictions["br_{PHI}_tautau_non-mass"].Clone(sfname.format(PHI=args.bsm_sm_like))
if args.bsm_sm_like in ['h', 'H']:
    bsm_predictions[sfname].Multiply(bsm_predictions["gsq_{PHI}_VV"])

# Ratio of the two types of scale factors for ggPHI and qqPHI
sfname = "r_sf_gg_{PHI}"
bsm_predictions[sfname] = bsm_predictions["sf_gg_{PHI}"].Clone(sfname.format(PHI=args.bsm_sm_like))
bsm_predictions[sfname].Divide(bsm_predictions["sf_gg_{PHI}_mass-corr"])

sfname = "r_sf_bb_{PHI}"
bsm_predictions[sfname] = bsm_predictions["sf_bb_{PHI}"].Clone(sfname.format(PHI=args.bsm_sm_like))
bsm_predictions[sfname].Divide(bsm_predictions["sf_bb_{PHI}_mass-corr"])

sfname = "r_sf_qq_{PHI}"
bsm_predictions[sfname] = bsm_predictions["sf_qq_{PHI}"].Clone(sfname.format(PHI=args.bsm_sm_like))
bsm_predictions[sfname].Divide(bsm_predictions["sf_qq_{PHI}_mass-corr"])

# Compute the contours for invalid mass values of SM-like Higgs boson
mh122_contours = plot.contourFromTH2(bsm_predictions["m_{PHI}_inverted"], (1-1./mass_borders[0]), 5, frameValue=1)
mh128_contours = plot.contourFromTH2(bsm_predictions["m_{PHI}"], mass_borders[1], 5, frameValue=1)

for graph in mh122_contours:
    if graph.GetN() > 5:
        graph.SetLineColor(R.kRed)
        graph.SetLineWidth(3)
        graph.SetFillColor(R.kRed)
        graph.SetFillStyle(3004)
        contour_graphs.setdefault("m_{PHI}_border", []).append((mass_borders[0], graph.Clone()))

for graph in mh128_contours:
    if graph.GetN() > 5:
        graph.SetLineColor(R.kRed)
        graph.SetLineWidth(3)
        graph.SetFillColor(R.kRed)
        graph.SetFillStyle(3004)
        contour_graphs.setdefault("m_{PHI}_border", []).append((mass_borders[1], graph.Clone()))

legend_mphi = R.TLegend(0.08,0.95,0.6,0.99)
legend_mphi.SetFillStyle(0)
legend_mphi.SetTextSize(0.03)
legend_mphi.AddEntry(contour_graphs["m_{PHI}_border"][0][1],"m_{PHI} #notin [122,128] GeV".replace("PHI","h_{1}" if args.bsm_sm_like == "H1" else args.bsm_sm_like),"F")

C.Clear()

out = R.TFile.Open(bsm_name + "_debug.root", "recreate")

for bsm_pred in bsm_predictions.values():
    print(bsm_pred)
    bsm_pred.Write()

out.Close()

# Compute the contours for required quantities 

contour_quantities = []

for key in bsm_predictions.keys():
    if "sf_" in key or "mass-only" in key or "non-mass" in key or key == "m_{PHI}" or "gsq_" in key:
        contour_quantities.append(key)

for key in contour_quantities:
    bsm_pred = bsm_predictions[key]
    contours = np.array(list(quantity_settings[key]["contours"].keys()))
    contour_graphs.setdefault(key, [])

    for cval in contours:
        contour_hist = bsm_predictions[key].Clone("conthist")
        contour_hist.SetContour(1, np.array([cval]))
        contour_hist.Draw('cont z list')
        C.Update()
        conts = R.gROOT.GetListOfSpecials().FindObject('contours')
        for cont in conts:
            for graph in cont:
                if graph.GetN() > 30:
                    graph.SetLineWidth(3)
                    contour_graphs[key].append((cval, graph.Clone()))
        C.Clear()

    # Resetting values to range borders, if magnitude too big
    for i_X in range(1, NXBins+1):
        for i_Y in range(1,NYBins+1):
            value = bsm_pred.GetBinContent(i_X, i_Y)
            if value <= quantity_settings[key]["range"][0]:
                bsm_pred.SetBinContent(i_X, i_Y, quantity_settings[key]["range"][0] + shift)
            if value >= quantity_settings[key]["range"][1]:
                bsm_pred.SetBinContent(i_X, i_Y, quantity_settings[key]["range"][1] - shift)

# Prepare plotting
haxis = bsm_predictions["m_{PHI}"].Clone("axis")
xtitle = "m_{A} [GeV]" if args.bsm_sm_like == "h" else "m_{H^{+}} [GeV]"
haxis.GetYaxis().SetTitleOffset(0.95)

latex = R.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)

if not os.path.isdir(os.path.join(args.plots, bsm_name)):
    os.makedirs(os.path.join(args.plots, bsm_name))


for key in contour_quantities:
        C.Clear()
        contour_legend = R.TLegend(0.6, 0.6, 0.9, 0.9)
        contour_legend.SetFillStyle(0)
        contour_legend.SetTextSize(0.04)
        haxis.SetMinimum(quantity_settings[key]["range"][0])
        haxis.SetMaximum(quantity_settings[key]["range"][1])
        haxis.SetTitle(";".join(["",xtitle,"tan#beta"]))
        if "EFT" in bsm_name:
            haxis.GetXaxis().SetRangeUser(91.5, haxis.GetXaxis().GetXmax())
        haxis.Draw("axis")
        bsm_predictions[key].GetZaxis().SetTitle(quantity_settings[key]["name"])
        bsm_predictions[key].GetZaxis().SetTitleOffset(1.4)
        bsm_predictions[key].Draw("colz same")
        current_level = None
        for level,graph in contour_graphs[key]:
                graph.SetLineColor(quantity_settings[key]["contours"][level])
                if current_level != level:
                    current_level = level
                    contour_legend.AddEntry(graph, str(current_level), "l")
                graph.Draw("C same")

        if key != "m_{PHI}":
            for level,graph in contour_graphs["m_{PHI}_border"]:
                    graph.Draw("C same")
                    graph.Draw("F same")

        contour_legend.Draw()
        if key != "m_{PHI}":
            legend_mphi.Draw()
        latex.DrawLatex(haxis.GetXaxis().GetXmax(), haxis.GetYaxis().GetXmax()+(haxis.GetYaxis().GetXmax()-haxis.GetYaxis().GetXmin())*0.02, bsm_model_names[bsm_name])
        C.Update()
        C.RedrawAxis()
        plotname = os.path.join(args.plots, bsm_name, key.format(PHI=args.bsm_sm_like))
        C.SaveAs(plotname + ".pdf")
        C.SaveAs(plotname + ".png")
