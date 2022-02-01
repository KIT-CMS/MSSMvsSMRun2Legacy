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

args = parser.parse_args()


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
    "m_A" : "m_{A} [GeV]",
    "m_H" : "m_{H} [GeV]",
    "m_h" : "m_{h} [GeV]",
    "m_H1" : "m_{h_{1}} [GeV]",
    "m_H2" : "m_{h_{2}} [GeV]",
    "m_H3" : "m_{h_{3}} [GeV]",
    "xs_gg_A" : "#sigma(ggA) [pb]",
    "xs_gg_H" : "#sigma(ggH) [pb]",
    "xs_gg_h" : "#sigma(ggh) [pb]",
    "xs_gg_H1" : "#sigma(gg^{}h_{1}) [pb]",
    "xs_gg_H2" : "#sigma(gg^{}h_{2}) [pb]",
    "xs_gg_H3" : "#sigma(gg^{}h_{3}) [pb]",
    "br_A_tautau" : "BR(A#rightarrow#tau#tau)",
    "br_H_tautau" : "BR(H#rightarrow#tau#tau)",
    "br_h_tautau" : "BR(h#rightarrow#tau#tau)",
    "br_H1_tautau" : "BR(^{}h_{1}#rightarrow#tau#tau)",
    "br_H2_tautau" : "BR(^{}h_{2}#rightarrow#tau#tau)",
    "br_H3_tautau" : "BR(^{}h_{3}#rightarrow#tau#tau)",
}

hists = {}

if not os.path.exists(model):
    os.makedirs(model)

c.Clear()
hists[model] = {}

smlike = "h"
additional = ["H", "A"]
if model == "mHH125":
    smlike = "H"
    additional = ["h", "A"]
elif model == "mh1125_CPV":
    smlike = "H1"
    additional = ["H2", "H3"]

mphi_parameter = "m_{H^{#pm}} [GeV]" if model in ["mHH125", "mh1125_CPV"] else "m_{A} [GeV]"

hists[model]["m_{smlike}".format(smlike=smlike)] = modelfile.Get("m_{smlike}".format(smlike=smlike))
hists[model]["m_{smlike}".format(smlike=smlike)].GetXaxis().SetTitle(mphi_parameter)
hists[model]["m_{smlike}".format(smlike=smlike)].GetYaxis().SetTitle("tan#beta")
hists[model]["m_{smlike}".format(smlike=smlike)].GetYaxis().SetTitleOffset(0.9)
hists[model]["m_{smlike}".format(smlike=smlike)].GetZaxis().SetTitleOffset(1.4)
hists[model]["m_{smlike}".format(smlike=smlike)].GetZaxis().SetTitle(quantity_dict["m_{smlike}".format(smlike=smlike)])
hists[model]["m_{smlike}".format(smlike=smlike)].SetMaximum(min(135., hists[model]["m_{smlike}".format(smlike=smlike)].GetMaximum()))
hists[model]["m_{smlike}".format(smlike=smlike)].SetMinimum(max(115., hists[model]["m_{smlike}".format(smlike=smlike)].GetMinimum()))
hists[model]["m_{smlike}".format(smlike=smlike)].SetContour(40)
hists[model]["m_{smlike}".format(smlike=smlike)].Draw("colz")
c.RedrawAxis()
c.Update()
text = r.TLatex()
text.SetTextAlign(11)
text.SetTextFont(42)
text.SetTextSize(0.05)
text.DrawLatex(hists[model]["m_{smlike}".format(smlike=smlike)].GetXaxis().GetXmin(), hists[model]["m_{smlike}".format(smlike=smlike)].GetYaxis().GetXmax() + (hists[model]["m_{smlike}".format(smlike=smlike)].GetYaxis().GetXmax() - hists[model]["m_{smlike}".format(smlike=smlike)].GetYaxis().GetXmin())*0.03, scenario_names[model])
plotname = "_".join([model, "m_{smlike}".format(smlike=smlike)])
c.SaveAs(os.path.join(model,plotname+".png"))
c.Print(os.path.join(model,plotname+".pdf"))

c.Clear()
hists[model]["xs_gg_{smlike}".format(smlike=smlike)] = modelfile.Get("xs_gg_{smlike}".format(smlike=smlike))
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetXaxis().SetTitle(mphi_parameter)
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().SetTitle("tan#beta")
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().SetTitleOffset(0.9)
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetZaxis().SetTitleOffset(1.4)
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetZaxis().SetTitle(quantity_dict["xs_gg_{smlike}".format(smlike=smlike)])
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].SetContour(40)
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].Draw("colz")
c.RedrawAxis()
c.Update()
text.DrawLatex(hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetXaxis().GetXmin(), hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().GetXmax() + (hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().GetXmax() - hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().GetXmin())*0.03, scenario_names[model])
plotname = "_".join([model, "xs_gg_{smlike}".format(smlike=smlike)])
c.SaveAs(os.path.join(model,plotname+".png"))
c.Print(os.path.join(model,plotname+".pdf"))
c.Clear()
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].SetMaximum(min(100., hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetMaximum()))
hists[model]["xs_gg_{smlike}".format(smlike=smlike)].Draw("colz")
c.RedrawAxis()
c.Update()
text.DrawLatex(hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetXaxis().GetXmin(), hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().GetXmax() + (hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().GetXmax() - hists[model]["xs_gg_{smlike}".format(smlike=smlike)].GetYaxis().GetXmin())*0.03, scenario_names[model])
c.SaveAs(os.path.join(model,plotname+"_restricted.png"))
c.Print(os.path.join(model,plotname+"_restricted.pdf"))

c.Clear()
hists[model]["br_{smlike}_tautau".format(smlike=smlike)] = modelfile.Get("br_{smlike}_tautau".format(smlike=smlike))
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetXaxis().SetTitle(mphi_parameter)
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().SetTitle("tan#beta")
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().SetTitleOffset(0.9)
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetZaxis().SetTitleOffset(1.4)
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetZaxis().SetTitle(quantity_dict["br_{smlike}_tautau".format(smlike=smlike)])
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].SetContour(40)
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].Draw("colz")
c.RedrawAxis()
c.Update()
text.DrawLatex(hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetXaxis().GetXmin(), hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().GetXmax() + (hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().GetXmax() - hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().GetXmin())*0.03, scenario_names[model])
plotname = "_".join([model, "br_{smlike}_tautau".format(smlike=smlike)])
c.SaveAs(os.path.join(model,plotname+".png"))
c.Print(os.path.join(model,plotname+".pdf"))
c.Clear()
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].SetMaximum(min(0.1, hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetMaximum()))
hists[model]["br_{smlike}_tautau".format(smlike=smlike)].Draw("colz")
c.RedrawAxis()
c.Update()
text.DrawLatex(hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetXaxis().GetXmin(), hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().GetXmax() + (hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().GetXmax() - hists[model]["br_{smlike}_tautau".format(smlike=smlike)].GetYaxis().GetXmin())*0.03, scenario_names[model])
c.SaveAs(os.path.join(model,plotname+"_restricted.png"))
c.Print(os.path.join(model,plotname+"_restricted.pdf"))
