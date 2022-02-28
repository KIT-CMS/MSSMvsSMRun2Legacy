#!/usr/bin/env python
import ROOT as r
from CombineHarvester.CombineTools.plotting import ModTDRStyle
from array import array
import argparse
import numpy as np
import os
import sys

r.gROOT.SetBatch()
ModTDRStyle()

parser = argparse.ArgumentParser(description="Script to plot quantities from MSSM benchmark ROOT files")
parser.add_argument("--model", help="Path to the model ROOT file to be studied")
parser.add_argument("--restricted", action="store_true", help="Flag to restrict ratio to [0.6, 1.4]")

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

bsm_higgs_names = {
    "H" : "H",
    "h" : "h",
    "H1" : "^{}h_{1}",
}

modelfile = r.TFile.Open(args.model, "read")
model = None
for name in scenario_names:
    if name + "_13.root" in os.path.basename(args.model):
        model = name
        break
if not model:
    print("Cannot determine model from file name")
    sys.exit(1)

quantity_dict = {
    "xs_gg_H" : "#sigma(gg#phi)",
    "xs_gg_h" : "#sigma(gg#phi)",
    "xs_gg_H1" : "#sigma(gg#phi)",
    "xs_bb_H" : "#sigma(bb#phi)",
    "xs_bb_h" : "#sigma(bb#phi)",
    "xs_bb_H1" : "#sigma(bb#phi)",
    "xs_vbf_H" : "#sigma(VBF#phi)",
    "xs_vbf_h" : "#sigma(VBF#phi)",
    "xs_vbf_H1" : "#sigma(VBF#phi)",
    "xs_hs_ZH" : "#sigma(Z#phi)",
    "xs_hs_Zh" : "#sigma(Z#phi)",
    "xs_hs_ZH1" : "#sigma(Z#phi)",
    "xs_hs_WH" : "#sigma(W#phi)",
    "xs_hs_Wh" : "#sigma(W#phi)",
    "xs_hs_WH1" : "#sigma(W#phi)",
    "xs_tth_H" : "#sigma(tt#phi)",
    "xs_tth_h" : "#sigma(tt#phi)",
    "xs_tth_H1" : "#sigma(tt#phi)",
    "br_H_tautau" : "BR(#phi#rightarrow#tau#tau)",
    "br_h_tautau" : "BR(#phi#rightarrow#tau#tau)",
    "br_H1_tautau" : "BR(#phi#rightarrow#tau#tau)",
    "m_H" : "m_{#phi}",
    "m_h" : "m_{#phi}",
    "m_H1" :"m_{#phi}",
}

hists = {}

if not os.path.exists(model):
    os.makedirs(model)

c.Clear()
hists[model] = {}

smlike = "h"
if model == "mHH125":
    smlike = "H"
elif model == "mh1125_CPV":
    smlike = "H1"

mphi_parameter = "m_{H^{#pm}} [GeV]" if model in ["mHH125", "mh1125_CPV"] else "m_{A} [GeV]"
nbinsx = modelfile.Get("m_"+smlike).GetXaxis().GetNbins()
nbinsy = modelfile.Get("m_"+smlike).GetYaxis().GetNbins()

zmin = 0.6
zmax = 1.4

quantities = ["xs_gg_{smlike}", "xs_bb_{smlike}", "xs_vbf_{smlike}", "xs_hs_Z{smlike}", "xs_hs_W{smlike}", "xs_tth_{smlike}", "br_{smlike}_tautau", "m_{smlike}"]

for q in quantities:
    histname = q.format(smlike="phi")+"_ratio"
    hists[histname] = modelfile.Get(q.format(smlike=smlike))
    if q == "m_{smlike}":
        hists[histname].Scale(1./125.38)
    else:
        hists[histname].Divide(modelfile.Get(q.format(smlike="HSM")))
    hists[histname].GetXaxis().SetTitle(mphi_parameter)
    hists[histname].GetYaxis().SetTitle("tan#beta")
    hists[histname].GetYaxis().SetTitleOffset(0.9)
    hists[histname].GetZaxis().SetTitleOffset(1.4)
    if q == "m_{smlike}":
        hists[histname].GetZaxis().SetTitle("m_{{smlike}}/(125.38 GeV)".replace("{smlike}",bsm_higgs_names[smlike]))
    else:
        hists[histname].GetZaxis().SetTitle("{smlike}/^{}H_{SM} ratio for ".replace("{smlike}",bsm_higgs_names[smlike])+quantity_dict[q.format(smlike=smlike)])
    if args.restricted:
        for i in range(nbinsx):
            for j in range(nbinsy):
                qval = hists[histname].GetBinContent(i+1,j+1)
                qval = qval if qval < zmax else zmax
                qval = qval if qval > zmin else zmin
                if qval in [zmin, zmax]:
                    hists[histname].SetBinContent(i+1,j+1, qval)
        hists[histname].SetMaximum(zmax)
        hists[histname].SetMinimum(zmin)
    hists[histname].SetContour(40)
    hists[histname].Draw("colz")
    c.RedrawAxis()
    c.Update()
    text = r.TLatex()
    text.SetTextAlign(11)
    text.SetTextFont(42)
    text.SetTextSize(0.05)
    text.DrawLatex(hists[histname].GetXaxis().GetXmin(), hists[histname].GetYaxis().GetXmax() + (hists[histname].GetYaxis().GetXmax() - hists[histname].GetYaxis().GetXmin())*0.03, scenario_names[model])
    if args.restricted:
        c.SaveAs(os.path.join(model,histname+"_BSMvsSM_restricted.png"))
        c.Print(os.path.join(model,histname+"_BSMvsSM_restricted.pdf"))
    else:
        c.SaveAs(os.path.join(model,histname+"_BSMvsSM.png"))
        c.Print(os.path.join(model,histname+"_BSMvsSM.pdf"))
