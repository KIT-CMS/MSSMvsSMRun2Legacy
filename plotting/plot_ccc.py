#!/usr/bin/env python
# The script is inspired by the macro https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/102x/test/plotting/cccPlot.cxx

import argparse
import os

import ROOT
ROOT.gROOT.SetBatch()
ROOT.PyConfig.IgnoreCommandLineOptions = True


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,
                        help="Input root file.")
    parser.add_argument("-p", "--poi",
                        type=str,
                        default="r_ggH",
                        help="POI that is fit.")
    parser.add_argument("-o", "--output",
                        default="ChannelCompatibilityCheck",
                        help="Name of the plot file")
    parser.add_argument("-r", "--frame-range",
                        type=lambda ranges: [float(edge.replace("m", "-")) for edge in ranges.split(",")],
                        default=[-10, 10],
                        help="Range of the fit values")
    parser.add_argument("--separate-nominal-fit",
                        type=str,
                        default=None,
                        help="Fallback solution, take nominal fit result from different fit")
    parser.add_argument("--mass",
                        type=str,
                        default=None,
                        help="Mass of the considered Higgs boson.")
    return parser.parse_args()


def main(args):
    # Read in input root file.
    if not os.path.exists(args.input):
        print("Input file {} does not exist.".format(args.input))
        raise ValueError
    infile = ROOT.TFile(args.input, "read")

    # Load the fit results
    if args.separate_nominal_fit is None:
        fit_nominal = infile.Get("fit_nominal")
    else:
        if not os.path.exists(args.separate_nominal_fit):
            print("Input file {} for separate fit does not exist.".format(args.input))
            raise ValueError
        infile_separate = ROOT.TFile(args.separate_nominal_fit, "read")
        fit_nominal = infile_separate.Get("fit_nominal")
    fit_alternate = infile.Get("fit_alternate")

    # Get POI result RooFit object
    res_poi = fit_nominal.floatParsFinal().find(args.poi)

    # Get number of channels from fit result
    prefix = "_ChannelCompatibilityCheck_{}_".format(args.poi)
    parameters = fit_alternate.floatParsFinal()
    num_pars = parameters.getSize()
    channels = [parameters.at(x).GetName() for x in range(num_pars) if parameters.at(x).GetName().startswith(prefix)]
    num_chans = len(channels)
    print(channels)

    # canv = ROOT.TCanvas("canv", "canv", 700, 1400)
    canv = ROOT.TCanvas("canv")
    canv.SetLeftMargin(0.2)
    canv.SetBottomMargin(0.13)
    canv.SetGridx(1)
    # Create the plot object.
    frame = ROOT.TH2F("frame",
                      ";best fit {};".format(args.poi.replace("_ggH", "_{ggH}")),
                      1,
                      # res_poi.getVal() + 5*res_poi.getAsymErrorLo(),
                      args.frame_range[0],
                      args.frame_range[1],
                      num_chans,
                      0,
                      num_chans
                      )

    print("Nominal fit: {} {}/+{}".format(res_poi.getVal(), res_poi.getAsymErrorLo(), res_poi.getAsymErrorHi()))
    # Fill TGraphAsymmErrors with fit values.
    points = ROOT.TGraphAsymmErrors(num_chans)
    for i, ch_ in enumerate(reversed(sorted(channels))):
        ri = parameters.find(ch_)
        points.SetPoint(i, ri.getVal(), i+0.5)
        points.SetPointError(i, -ri.getAsymErrorLo(), ri.getAsymErrorHi(), 0, 0)
        print("Alternate fit: ", ch_, ri.getVal())
        frame.GetYaxis().SetBinLabel(i+1, ch_.replace(prefix, ""))

    points.SetLineColor(ROOT.kBlack)
    points.SetLineWidth(3)
    points.SetMarkerStyle(21)
    # frame.GetXaxis().SetNdivisions(505)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.06)
    frame.Draw()

    ROOT.gStyle.SetOptStat(0)
    globalFitBand = ROOT.TBox(res_poi.getVal()+res_poi.getAsymErrorLo(), 0, res_poi.getVal()+res_poi.getAsymErrorHi(), num_chans)
    globalFitBand.SetFillStyle(3013)
    globalFitBand.SetFillColor(65)
    globalFitBand.SetLineStyle(0)
    globalFitBand.Draw("same")

    globalFitLine = ROOT.TLine(res_poi.getVal(), 0, res_poi.getVal(), num_chans)
    globalFitLine.SetLineWidth(4)
    globalFitLine.SetLineColor(214)
    globalFitLine.Draw("same")
    points.Draw("P SAME")

    # legend = ROOT.TLegend(0.6, 0.13, 0.9, 0.33)
    legend = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
    # legend = ROOT.TLegend(0.2, 0.13, 0.5, 0.33)
    legend.AddEntry(globalFitLine, "Global Best Fit", "l")
    legend.AddEntry(globalFitBand, "Global Best Fit #pm 1 #sigma", "f")
    legend.Draw()

    limit_tree = infile.Get("limit")
    limit_tree.GetEntry(0)
    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(ROOT.kBlack)
    latex2.SetTextSize(0.04)
    begin_left = 0.245
    # latex2.DrawLatex(begin_left, 0.360, "#chi^{{2}}-like = {:.3f}".format(limit_tree.limit))
    latex2.DrawLatex(begin_left, 0.640, "#chi^{{2}}-like = {:.3f}".format(limit_tree.limit))

    # Draw mass of Higgs boson outside of frame
    if args.mass is not None:
        latex2 = ROOT.TLatex()
        latex2.SetNDC()
        latex2.SetTextAngle(0)
        latex2.SetTextColor(ROOT.kBlack)
        latex2.SetTextSize(0.04)
        begin_left = 0.205
        latex2.DrawLatex(begin_left, 0.920, "m_{#phi} = %s GeV" % args.mass)

    canv.Print(args.output + ".png", "png")
    canv.Print(args.output + ".pdf", "pdf")
    return


if __name__ == "__main__":
    args = parse_args()
    main(args)
