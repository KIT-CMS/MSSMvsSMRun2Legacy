#!/usr/bin/env python

import ROOT as r
from CombineHarvester.CombineTools.plotting import ModTDRStyle
from array import array
import numpy as np
import argparse
import os
import sys

r.gROOT.SetBatch()
ModTDRStyle()

parser = argparse.ArgumentParser(description="Script to plot differences for previous (on Twiki) and current (on zenodo) model releases for MSSM benchmarks")
parser.add_argument("--old-model", help="Path to the model ROOT file to be studied, for old release on Twiki")
parser.add_argument("--new-model", help="Path to the model ROOT file to be studied, for new release on zenodo")

args = parser.parse_args()


r.TColor.CreateGradientColorTable(3,
                                  array("d",[0.  ,0.5,  1.]),
                                  array("d",[0.  , 1.,  1.]),
                                  array("d",[0.35, 1.,0.65]),
                                  array("d",[1.  , 1.,  0.]),
                                  10000,  1.0)

c = r.TCanvas("c", "c", 1600, 1600)
c.SetLeftMargin(1.3)
c.SetRightMargin(0.2)
c.SetTopMargin(1.1)
c.cd()

scenario_names = {
    "mh125" : "M_{h}^{125}",
    "mh125EFT" : "M_{h,EFT}^{125}",
    "mh125_lc" : "M_{h}^{125}(#tilde{#chi})",
    "mh125EFT_lc" : "M_{h,EFT}^{125}(#tilde{#chi})",
    "mh125_ls" : "M_{h}^{125}(#tilde{#tau})",
    "mh125_align" : "M_{h}^{125}(alignment)",
    "mHH125" : "M_{H}^{125}",
    "mh1125_CPV" : "M_{h_{1}}^{125}(CPV)",
    "mh125_muneg_1" : "M_{h}^{125 ^{}#mu_{1}#minus}",
    "mh125_muneg_2" : "M_{h}^{125 ^{}#mu_{2}#minus}",
    "mh125_muneg_3" : "M_{h}^{125 ^{}#mu_{3}#minus}",
    "hMSSM" : "hMSSM",
}

files = {
    "previous" : r.TFile.Open(args.old_model, "read"),
    "current" : r.TFile.Open(args.new_model, "read"),
}

model = None
old_model = None
new_model = None
for name in scenario_names:
    if name + "_13_old.root" in os.path.basename(args.old_model):
        old_model = name
        break
for name in scenario_names:
    if name + "_13.root" in os.path.basename(args.new_model):
        new_model = name
        break
if not old_model == new_model:
    print("Two different models where chosen")
    sys.exit(1)
else:
    model = new_model
    
if not model:
    print("Cannot determine model from file name")
    sys.exit(1)


smlike = "h"
additional = ["H", "A"]
if model == "mHH125":
    smlike = "H"
    additional = ["h", "A"]
elif model == "mh1125_CPV":
    smlike = "H1"
    additional = ["H2", "H3"]

mphi_parameter = "m_{H^{#pm}} [GeV]" if model in ["mHH125", "mh1125_CPV"] else "m_{A} [GeV]"

zmin = 0.6
zmax = 1.4

nbinsx = files["previous"].Get("m_"+smlike).GetXaxis().GetNbins()
nbinsy = files["previous"].Get("m_"+smlike).GetYaxis().GetNbins()

quantity_dict = {
    "xs_gg_A" : "#sigma(ggA)",
    "xs_gg_H" : "#sigma(ggH)",
    "xs_gg_h" : "#sigma(ggh)",
    "xs_gg_H1" : "#sigma(gg^{}h_{1}) [pb]",
    "xs_gg_H2" : "#sigma(gg^{}h_{2}) [pb]",
    "xs_gg_H3" : "#sigma(gg^{}h_{3}) [pb]",
    "xs_gg_HSM" : "#sigma(gg^{}H_{SM})",
    "xs_bb_A" : "#sigma(bbA)",
    "xs_bb_H" : "#sigma(bbH)",
    "xs_bb_h" : "#sigma(bbh)",
    "xs_bb_H1" : "#sigma(bb^{}h_{1}) [pb]",
    "xs_bb_H2" : "#sigma(bb^{}h_{2}) [pb]",
    "xs_bb_H3" : "#sigma(bb^{}h_{3}) [pb]",
    "xs_bb_HSM" : "#sigma(bb^{}H_{SM})",
    "br_A_tautau" : "BR(A#rightarrow#tau#tau)",
    "br_H_tautau" : "BR(H#rightarrow#tau#tau)",
    "br_h_tautau" : "BR(h#rightarrow#tau#tau)",
    "br_H1_tautau" : "BR(^{}h_{1}#rightarrow#tau#tau)",
    "br_H2_tautau" : "BR(^{}h_{2}#rightarrow#tau#tau)",
    "br_H3_tautau" : "BR(^{}h_{3}#rightarrow#tau#tau)",
    "br_HSM_tautau" : "BR(^{}H_{SM}#rightarrow#tau#tau)",
    "br_A_bb" : "BR(A#rightarrowbb)",
    "br_H_bb" : "BR(H#rightarrowbb)",
    "br_h_bb" : "BR(h#rightarrowbb)",
    "br_H1_bb" : "BR(^{}h_{1}#rightarrowbb)",
    "br_H2_bb" : "BR(^{}h_{2}#rightarrowbb)",
    "br_H3_bb" : "BR(^{}h_{3}#rightarrowbb)",
    "br_HSM_bb" : "BR(^{}H_{SM}#rightarrowbb)",
    "rescale_gt_A" : "^{}Y_{t}(A)",
    "rescale_gt_H" : "^{}Y_{t}(H)",
    "rescale_gt_h" : "^{}Y_{t}(h)",
    "rescale_gt_H1" : "^{}Y_{t}(^{}h_{1})",
    "rescale_gt_H2" : "^{}Y_{t}(^{}h_{2})",
    "rescale_gt_H3" : "^{}Y_{t}(^{}h_{3})",
    "rescale_gb_A" : "^{}Y_{b}(A)",
    "rescale_gb_H" : "^{}Y_{b}(H)",
    "rescale_gb_h" : "^{}Y_{b}(h)",
    "rescale_gb_H1" : "^{}Y_{b}(^{}h_{1})",
    "rescale_gb_H2" : "^{}Y_{b}(^{}h_{2})",
    "rescale_gb_H3" : "^{}Y_{b}(^{}h_{3})",
    "width_A" : "^{}#Gamma_{tot}(A)",
    "width_H" : "^{}#Gamma_{tot}(H)",
    "width_h" : "^{}#Gamma_{tot}(h)",
    "width_H1" : "^{}#Gamma_{tot}(^{}h_{1})",
    "width_H2" : "^{}#Gamma_{tot}(^{}h_{2})",
    "width_H3" : "^{}#Gamma_{tot}(^{}h_{3})",
    "width_HSM" : "^{}#Gamma_{tot}(^{}H_{SM})",
    "int_gg_tautau_H2" : "^{}#eta(gg#rightarrow^{}h_{2}#rightarrow#tau#tau)",
    "int_gg_tautau_H3" : "^{}#eta(gg#rightarrow^{}h_{3}#rightarrow#tau#tau)",
    "int_bb_tautau_H2" : "^{}#eta(bb#rightarrow^{}h_{2}#rightarrow#tau#tau)",
    "int_bb_tautau_H3" : "^{}#eta(bb#rightarrow^{}h_{3}#rightarrow#tau#tau)",
}

BSM = {}
for phi in additional + [smlike]:
    BSM[phi] = [q.format(phi=phi) for q in ["xs_gg_{phi}", "xs_bb_{phi}", "rescale_gb_{phi}", "rescale_gt_{phi}", "br_{phi}_tautau", "br_{phi}_bb", "width_{phi}"]]
    if model == "mh1125_CPV" and phi in additional:
        BSM[phi].extend([q.format(phi=phi) for q in ["int_gg_tautau_{phi}", "int_bb_tautau_{phi}"]])

SM = [{"previous" : q + "_SM", "current" : q.replace("_"+smlike, "_HSM")} for q in BSM[smlike] if not "rescale_" in q and not "int_" in q]

if not os.path.exists(model):
    os.makedirs(model)

hists = {}
# Difference between BSM quantities
for higgs,quantities in BSM.items():
    hists[higgs] = {}
    for q in quantities:
        c.Clear()
        hists[higgs][q] = {}
        valid_sm_inputs = True
        for release,inputfile in files.items():
            hists[higgs][q][release] = inputfile.Get(q)
            if not hists[higgs][q][release]:
                print("Skipping "+q+" because quantity missing for release "+release+" of model "+model)
                valid_sm_inputs = False
        if not valid_sm_inputs:
            continue
        hists[higgs][q]["ratio"] = hists[higgs][q]["current"].Clone()
        hists[higgs][q]["ratio"].Divide(hists[higgs][q]["previous"])
        hists[higgs][q]["ratio"].Draw("colz")
        hists[higgs][q]["ratio"].GetXaxis().SetTitle(mphi_parameter)
        hists[higgs][q]["ratio"].GetYaxis().SetTitle("tan#beta")
        hists[higgs][q]["ratio"].GetYaxis().SetTitleOffset(0.9)
        hists[higgs][q]["ratio"].GetZaxis().SetTitle("current/previous for {quantity}".format(quantity=quantity_dict[q]))
        hists[higgs][q]["ratio"].GetZaxis().SetTitleOffset(1.1)
        for i in range(nbinsx):
            for j in range(nbinsy):
                qval = hists[higgs][q]["ratio"].GetBinContent(i+1,j+1)
                qval = qval if qval < zmax else zmax
                qval = qval if qval > zmin else zmin
                if qval in [zmin, zmax]:
                    hists[higgs][q]["ratio"].SetBinContent(i+1,j+1, qval)
        hists[higgs][q]["ratio"].SetMinimum(zmin)
        hists[higgs][q]["ratio"].SetMaximum(zmax)
        hists[higgs][q]["ratio"].SetContour(40)
        c.RedrawAxis()
        c.Update()
        text = r.TLatex()
        text.SetTextAlign(11)
        text.SetTextFont(42)
        text.SetTextSize(0.05)
        text.DrawLatex(hists[higgs][q]["ratio"].GetXaxis().GetXmin(), hists[higgs][q]["ratio"].GetYaxis().GetXmax() + (hists[higgs][q]["ratio"].GetYaxis().GetXmax() - hists[higgs][q]["ratio"].GetYaxis().GetXmin())*0.03, scenario_names[model])
        c.SaveAs(os.path.join(model,q+"_comparison.png"))
        c.Print(os.path.join(model,q+"_comparison.pdf"))


# Difference between SM quantities
hists["HSM"] = {}
for qSM in SM:
    q = qSM["current"]
    c.Clear()
    hists["HSM"][q] = {}
    valid_sm_inputs = True
    for release,inputfile in files.items():
        hists["HSM"][q][release] = inputfile.Get(qSM[release])
        if not hists["HSM"][q][release]:
            print("Skipping "+q+" because quantity missing for release "+release+" of model "+model)
            valid_sm_inputs = False
    if not valid_sm_inputs:
        continue
    hists["HSM"][q]["ratio"] = hists["HSM"][q]["current"].Clone()
    hists["HSM"][q]["ratio"].Divide(hists["HSM"][q]["previous"])
    hists["HSM"][q]["ratio"].GetXaxis().SetTitle(mphi_parameter)
    hists["HSM"][q]["ratio"].GetYaxis().SetTitle("tan#beta")
    hists["HSM"][q]["ratio"].GetYaxis().SetTitleOffset(0.9)
    hists["HSM"][q]["ratio"].GetZaxis().SetTitle("current/previous for {quantity}".format(quantity=quantity_dict[q]))
    hists["HSM"][q]["ratio"].GetZaxis().SetTitleOffset(1.1)
    for i in range(nbinsx):
        for j in range(nbinsy):
            qval = hists["HSM"][q]["ratio"].GetBinContent(i+1,j+1)
            qval = qval if qval < zmax else zmax
            qval = qval if qval > zmin else zmin
            if qval in [zmin, zmax]:
                hists["HSM"][q]["ratio"].SetBinContent(i+1,j+1, qval)
    hists["HSM"][q]["ratio"].SetMinimum(zmin)
    hists["HSM"][q]["ratio"].SetMaximum(zmax)
    hists["HSM"][q]["ratio"].SetContour(40)
    hists["HSM"][q]["ratio"].Draw("colz")
    c.RedrawAxis()
    c.Update()
    text = r.TLatex()
    text.SetTextAlign(11)
    text.SetTextFont(42)
    text.SetTextSize(0.05)
    text.DrawLatex(hists["HSM"][q]["ratio"].GetXaxis().GetXmin(), hists["HSM"][q]["ratio"].GetYaxis().GetXmax() + (hists["HSM"][q]["ratio"].GetYaxis().GetXmax() - hists["HSM"][q]["ratio"].GetYaxis().GetXmin())*0.03, scenario_names[model])
    c.SaveAs(os.path.join(model,q+"_comparison.png"))
    c.Print(os.path.join(model,q+"_comparison.pdf"))

if not model == "mh1125_CPV":
    # Difference between angle alpha
    c.Clear()
    hists["alpha"] = {}
    hists["alpha"]["current"] = files["current"].Get("alpha")
    hists["alpha"]["previous"] = files["previous"].Get("m_h").Clone()
    hists["rescale_gt_H"] = inputfile.Get("rescale_gt_H")
    for i in range(nbinsx):
        for j in range(nbinsy):
            tanb = hists["alpha"]["previous"].GetYaxis().GetBinLowEdge(j+1)
            yt_H = hists["rescale_gt_H"].GetBinContent(i+1, j+1)
            beta = np.arctan(tanb)
            hists["alpha"]["previous"].SetBinContent(i+1, j+1, np.arcsin(yt_H*np.sin(beta)))
    hists["alpha"]["ratio"] = hists["alpha"]["current"].Clone()
    hists["alpha"]["ratio"].Divide(hists["alpha"]["previous"])
    hists["alpha"]["ratio"].GetXaxis().SetTitle(mphi_parameter)
    hists["alpha"]["ratio"].GetYaxis().SetTitle("tan#beta")
    hists["alpha"]["ratio"].GetYaxis().SetTitleOffset(0.9)
    hists["alpha"]["ratio"].GetZaxis().SetTitle("current/previous for #alpha")
    hists["alpha"]["ratio"].GetZaxis().SetTitleOffset(1.1)
    for i in range(nbinsx):
        for j in range(nbinsy):
            qval = hists["alpha"]["ratio"].GetBinContent(i+1,j+1)
            qval = qval if qval < zmax else zmax
            qval = qval if qval > zmin else zmin
            if qval in [zmin, zmax]:
                hists["alpha"]["ratio"].SetBinContent(i+1,j+1, qval)
    hists["alpha"]["ratio"].SetMinimum(zmin)
    hists["alpha"]["ratio"].SetMaximum(zmax)
    hists["alpha"]["ratio"].SetContour(40)
    hists["alpha"]["ratio"].Draw("colz")
    c.RedrawAxis()
    c.Update()
    text = r.TLatex()
    text.SetTextAlign(11)
    text.SetTextFont(42)
    text.SetTextSize(0.05)
    text.DrawLatex(hists["alpha"]["ratio"].GetXaxis().GetXmin(), hists["alpha"]["ratio"].GetYaxis().GetXmax() + (hists["alpha"]["ratio"].GetYaxis().GetXmax() - hists["alpha"]["ratio"].GetYaxis().GetXmin())*0.03, scenario_names[model])
    c.SaveAs(os.path.join(model,"alpha_comparison.png"))
    c.Print(os.path.join(model,"alpha_comparison.pdf"))


    # Difference between VBF+VH cross-section SF
    c.Clear()
    hists["SF"] = {}
    hists["SF"]["previous"] = files["previous"].Get("m_h").Clone()
    for i in range(nbinsx):
        for j in range(nbinsy):
            tanb = hists["SF"]["previous"].GetYaxis().GetBinLowEdge(j+1)
            beta = np.arctan(tanb)
            alpha = hists["alpha"]["current"].GetBinContent(i+1, j+1)
            if model == "mHH125":
                hists["SF"]["previous"].SetBinContent(i+1,j+1,np.cos(beta - alpha)**2)
            else:
                hists["SF"]["previous"].SetBinContent(i+1,j+1,np.sin(beta - alpha)**2)

    hists["SF"]["current_numerator"] = files["current"].Get("xs_vbf_h")
    hists["SF"]["current_numerator"].Add(files["current"].Get("xs_hs_Wh"))
    hists["SF"]["current_numerator"].Add(files["current"].Get("xs_hs_Zh"))
    hists["SF"]["current_denominator"] = files["current"].Get("xs_vbf_HSM")
    hists["SF"]["current_denominator"].Add(files["current"].Get("xs_hs_WHSM"))
    hists["SF"]["current_denominator"].Add(files["current"].Get("xs_hs_ZHSM"))
    hists["SF"]["current"] = hists["SF"]["current_numerator"].Clone()
    hists["SF"]["current"].Divide(hists["SF"]["current_denominator"])
    hists["SF"]["ratio"] = hists["SF"]["current"].Clone()
    hists["SF"]["ratio"].Divide(hists["SF"]["previous"])
    hists["SF"]["ratio"].GetXaxis().SetTitle("m_{A} [GeV]")
    hists["SF"]["ratio"].GetYaxis().SetTitle("tan#beta")
    hists["SF"]["ratio"].GetYaxis().SetTitleOffset(0.9)
    if model == "mHH125":
        hists["SF"]["ratio"].GetZaxis().SetTitle("#sigma((VBF#plusV)H)/#sigma((VBF#plusV)^{}H_{SM})/cos^{2}(#beta#minus#alpha)")
    else:
        hists["SF"]["ratio"].GetZaxis().SetTitle("#sigma((VBF#plusV)h)/#sigma((VBF#plusV)^{}H_{SM})/sin^{2}(#beta#minus#alpha)")
    hists["SF"]["ratio"].GetZaxis().SetTitleOffset(1.1)
    for i in range(nbinsx):
        for j in range(nbinsy):
            qval = hists["SF"]["ratio"].GetBinContent(i+1,j+1)
            qval = qval if qval < zmax else zmax
            qval = qval if qval > zmin else zmin
            if qval in [zmin, zmax]:
                hists["SF"]["ratio"].SetBinContent(i+1,j+1, qval)
    hists["SF"]["ratio"].SetMinimum(zmin)
    hists["SF"]["ratio"].SetMaximum(zmax)
    hists["SF"]["ratio"].SetContour(40)
    hists["SF"]["ratio"].Draw("colz")
    c.RedrawAxis()
    c.Update()
    text = r.TLatex()
    text.SetTextAlign(11)
    text.SetTextFont(42)
    text.SetTextSize(0.05)
    text.DrawLatex(hists["SF"]["ratio"].GetXaxis().GetXmin(), hists["SF"]["ratio"].GetYaxis().GetXmax() + (hists["SF"]["ratio"].GetYaxis().GetXmax() - hists["SF"]["ratio"].GetYaxis().GetXmin())*0.03, scenario_names[model])
    c.SaveAs(os.path.join(model,"SF_comparison.png"))
    c.Print(os.path.join(model,"SF_comparison.pdf"))
