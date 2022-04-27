#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import Dumbledraw.dumbledraw as dd
#import Dumbledraw.rootfile_parser_inputshapes as rootfile_parser
import Dumbledraw.rootfile_parser as rootfile_parser
import Dumbledraw.styles as styles

import argparse
import copy
import yaml
import os
from array import array

import logging
logger = logging.getLogger("")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=
        "Plot categories using Dumbledraw from shapes produced by shape-producer module."
    )
    parser.add_argument(
        "-l", "--linear", action="store_true", help="Enable linear x-axis")
    parser.add_argument(
        "-c",
        "--channels",
        nargs="+",
        type=str,
        required=True,
        help="Channels")
    parser.add_argument("-e", "--era", type=str, required=True, help="Era")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="ROOT file with shapes of processes")
    parser.add_argument(
        "--control-variable",
        type=str,
        default=None,
        help="Enable plotting goodness of fit shapes for given variable")
    parser.add_argument(
        "--png", action="store_true", help="Save plots in png format")
    parser.add_argument(
        "--normalize-by-bin-width",
        action="store_true",
        help="Normelize plots by bin width")
    parser.add_argument(
        "--fake-factor",
        action="store_true",
        help="Fake factor estimation method used")
    parser.add_argument(
        "--embedding",
        action="store_true",
        help="Fake factor estimation method used")
    parser.add_argument(
        "-o", "--output-dir",
        help="Output directory for the plots.")
    parser.add_argument(
        "-m", "--mass",
        type=str,
        default=None,
        help="Higgs boson mass displayed in the legend.")
    parser.add_argument(
        "--cross-section-ggh", "--x-sec-ggh",
        default=None,
        type=str,
        help="Cross sections displayed in the legend.")
    parser.add_argument(
        "--cross-section-bbh", "--x-sec-bbh",
        default=None,
        type=str,
        help="Cross sections displayed in the legend.")
    parser.add_argument(
        "--tanbeta",
        default=None,
        type=str,
        help="Higgs boson mass displayed on the legend.")
    parser.add_argument(
        "--control-region",
        action="store_true",
        help="Skip signal categories")
    parser.add_argument(
        "--model-independent",
        action="store_true",
        help="Plot shapes from model independent analysis.")
    parser.add_argument(
        "--blinded",
        action="store_true",
        help="Do not draw data.")
    parser.add_argument(
        "--x-range",
        type=lambda xranges: [float(edge) for edge in xranges.split(',')],
        default=None,
        help="Smaller x-range used in the plot to zoom into problematic regions")
    parser.add_argument(
        "--combine-backgrounds",
        action="store_true",
        help="Combine minor backgrounds to single shape")
    parser.add_argument(
        "--b-only-fit",
        type=str,
        default=None,
        help="Add total background shape of b-only fit to ratio plot."
    )
    return parser.parse_args()


def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


def rebin_hist_for_logX(hist, xlow=1.0):
    bins = []
    contents= []
    errors = []
    # Get low bin edges, contents and errors of all bins.
    for i in range(1, hist.GetNbinsX()+1):
        bins.append(hist.GetBinLowEdge(i))
        contents.append(hist.GetBinContent(i))
        errors.append(hist.GetBinError(i))
    # Add low edge of overflow bin as high edge of last regular bin.
    bins.append(hist.GetBinLowEdge((hist.GetNbinsX()+1)))
    # Check if first bin extents to zero and if it does change it to larger value
    if bins[0] == 0.:
        bins[0] = xlow
    # Create the new histogram
    bin_arr = array("f", bins)
    hist_capped = ROOT.TH1F(hist.GetName(), hist.GetTitle(), len(bin_arr)-1, bin_arr)
    for i in range(0, hist_capped.GetNbinsX()):
        hist_capped.SetBinContent(i+1, contents[i])
        hist_capped.SetBinError(i+1, errors[i])
    return hist_capped


def main(args):
    if args.control_variable is None:
        channel_categories = {
            #"et": ["nobtag_tight", "btag_tight", "nobtag_loosemt", "nobtag_tight"]
            "et": ["10", "11", "12", "13", "14", "15", "16", "17", "18", "32", "33",  "35", "36"],
            #"mt": ["nobtag_tight", "btag_tight", "nobtag_loosemt", "nobtag_tight"]
            "mt": ["10", "11", "12", "13", "14", "15", "16", "17", "18", "32", "33",  "35", "36"],
            #"tt": ["nobtag", "btag"]
            "tt": ["10", "11", "12", "13", "14", "15", "16", "17", "32", "35"],
            "em": ["2", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "32", "33", "34", "35", "36", "37"]
        }
        if args.model_independent:
            channel_categories = {
                #"et": ["nobtag_tight", "btag_tight", "nobtag_loosemt", "nobtag_tight"]
                "et": ["32", "33",  "35", "36"],
                #"mt": ["nobtag_tight", "btag_tight", "nobtag_loosemt", "nobtag_tight"]
                "mt": ["32", "33",  "35", "36"],
                #"tt": ["nobtag", "btag"]
                "tt": ["32", "35"],
                "em": ["2","32", "33", "34", "35", "36", "37"]
            }
            channel_categories = {
                "et": ["32", "35"],
                "mt": ["32", "35"],
                "tt": ["32", "35"],
                # "em": ["32", "33", "35", "36"],
                "em": ["33", "36"],
                "lt": ["32", "35"],
            }
    else:
        channel_categories = {
            "et": ["100"],
            "mt": ["100"],
            "tt": ["100"],
            "em": ["100"]
        }
    channel_dict = {
        "ee": "ee",
        "em": "e#mu",
        "et": "e#tau_{h}",
        "mm": "#mu#mu",
        "mt": "#mu#tau_{h}",
        "lt": "#mu_{}#tau_{h}+e_{}#tau_{h}",
        # "lt": "#ell#tau_{h}",
        "tt": "#tau_{h}#tau_{h}"
    }
    if args.control_variable != None:
        category_dict = {"100": "inclusive"}
    else:
        category_dict = {
            "et": {
                "1": "2D ggH/qqH category",
                "13": "tt",
                "15": "zll",
                "16": "misc",
                "20": "Genuine #tau",
                "21": "Jet #rightarrow #tau_{h}",
                "32": "No b tag Tight-m_{T}",
                "33": "No b tag Loose-m_{T}",
                "35": "b tag Tight-m_{T}",
                "36": "b tag Loose-m_{T}",
                },
            "mt": {
                "1": "2D ggH/qqH category",
                "13": "tt",
                "15": "zll",
                "16": "misc",
                "20": "Genuine #tau",
                "21": "Jet #rightarrow #tau_{h}",
                "32": "No b tag Tight-m_{T}",
                "33": "No b tag Loose-m_{T}",
                "35": "b tag Tight-m_{T}",
                "36": "b tag Loose-m_{T}",
                },
            "lt": {
                "1": "2D ggH/qqH category",
                "13": "tt",
                "15": "zll",
                "16": "misc",
                "20": "Genuine #tau",
                "21": "Jet #rightarrow #tau_{h}",
                "32": "No b tag Tight-m_{T}",
                "33": "No b tag Loose-m_{T}",
                "35": "b tag Tight-m_{T}",
                "36": "b tag Loose-m_{T}",
                },
            "tt": {
                "1": "2D ggH/qqH category",
                "16": "misc",
                "20": "Genuine #tau",
                "21": "Jet #rightarrow #tau_{h}",
                "32": "No b tag",
                "35": "b tag",
                },
            "em": {
                "1": "2D ggH/qqH category",
                "2": "D_{#lower[-0.2]{#zeta}} < -35 GeV",
                "13": "tt",
                "14": "qcd",
                "16": "misc",
                "19": "diboson",
                "20": "Genuine #tau",
                "32": "No b tag, High-D_{#lower[-0.2]{#zeta}}",
                "33": "No b tag, Medium-D_{#lower[-0.2]{#zeta}}",
                "34": "No b tag, Low-D_{#lower[-0.2]{#zeta}}",
                "35": "b tag, High-D_{#lower[-0.2]{#zeta}}",
                "36": "b tag, Medium-D_{#lower[-0.2]{#zeta}}",
                "37": "b tag, Low-D_{#lower[-0.2]{#zeta}}",
                },

        }
    if args.linear == True:
        split_value = 0
    else:
        if args.normalize_by_bin_width:
            split_value = 1
        else:
            split_value = 101

    split_dict = {c: split_value for c in ["et", "mt", "tt", "em", "lt"]}

    bkg_processes = [
        "HSM", "VVL", "TTL", "ZL", "jetFakes", "EMB"
    ]
    if args.combine_backgrounds:
        bkg_processes = [
            "other", "TTL", "jetFakes", "EMB"
        ]
    if not args.fake_factor and args.embedding:
        bkg_processes = [
            "QCD", "VVJ", "VVL", "W", "TTJ", "TTL", "ZJ", "ZL", "EMB"
        ]
    if not args.embedding and args.fake_factor:
        bkg_processes = [
            "VVT", "VVJ", "TTT", "TTJ", "ZJ", "ZL", "jetFakes", "ZTT"
        ]
    if not args.embedding and not args.fake_factor:
        bkg_processes = [
            "QCD", "W", "VVJ", "VVL", "VVT", "TTJ", "TTL", "TTT", "ZJ", "ZL", "ZTT"
#            "QCD", "VVT", "VVJ", "W", "TTT", "TTJ", "ZJ", "ZL", "ZTT"
        ]
    all_bkg_processes = [b for b in bkg_processes]
    legend_bkg_processes = copy.deepcopy(bkg_processes)
    legend_bkg_processes.reverse()

    if "2016" in args.era:
        era = "_2016"
    elif "2017" in args.era:
        era = "_2017"
    elif "2018" in args.era:
        era = "_2018"
    elif "combined" in args.era:
        era = ""
    else:
        logger.critical("Era {} is not implemented.".format(args.era))
        raise Exception

    plots = []
    for channel in args.channels:
        if "em" in channel:
            if not args.embedding:
                bkg_processes = [
                    "HSM", "QCDMC", "VVT", "VVL", "W", "TTT", "TTL", "ZL", "ZTT"
                ]
            else:
                bkg_processes = [
                    "HSM", "QCD", "EWK", "TTL", "ZL", "EMB"
                ]
                if args.combine_backgrounds:
                    bkg_processes = [
                        "other", "QCD", "TTL", "EMB"
                    ]

        for category in channel_categories[channel]:
            # Set split value to 101 for no b-tag categories to make
            # transition between scales smoother
            if int(category) < 35 and int(category) > 10:
                split_dict[channel] += 100
            if args.control_region and category != "2":
                continue
            rootfile = rootfile_parser.Rootfile_parser(args.input, mode="CombineHarvesterMerged")
            legend_bkg_processes = copy.deepcopy(bkg_processes)
            legend_bkg_processes.reverse()
            # create plot
            if args.linear == True:
                plot = dd.Plot(
                    [0.3, [0.3, 0.28]], "ModTDR", r=0.04, l=0.14, width=600)
            else:
                plot = dd.Plot(
                    [0.55, [0.3, 0.28]], style="ModTDR", invert_pad_creation=True, r=0.04, l=0.14, width=600)

            # get background histograms
            for process in bkg_processes:
                if process in ["jetFakes", "jetFakesEMB"] and channel == "tt":
                    jetfakes_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, process).Clone(), xlow=30.)
                    jetfakes_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "wFakes"), xlow=30.))
                    plot.add_hist(jetfakes_hist, process, "bkg")
                elif process in ["HSM"]:
                    if channel == "em":
                        hsm_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "ggH125").Clone(), xlow=30.)
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "qqH125"), xlow=30.))
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "bbH125"), xlow=30.))
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ggHWW125"), xlow=30.))
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "qqHWW125"), xlow=30.))
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ZHWW125"), xlow=30.))
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "WHWW125"), xlow=30.))
                        plot.add_hist(hsm_hist, process, "bkg")
                    else:
                        hsm_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "ggH125").Clone(), xlow=30.)
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "qqH125"), xlow=30.))
                        hsm_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "bbH125"), xlow=30.))
                        plot.add_hist(hsm_hist, process, "bkg")
                elif process == "EWK" and channel == "em":
                    ewk_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "W").Clone(), xlow=30.)
                    ewk_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "VVL"), xlow=30.))
                    plot.add_hist(ewk_hist, process, "bkg")
                elif process in ["HWW"]:
                    hww_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "ggHWW125").Clone(), xlow=30.)
                    hww_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "qqHWW125"), xlow=30.))
                    hww_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ZHWW125"), xlow=30.))
                    hww_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "WHWW125"), xlow=30.))
                    plot.add_hist(hww_hist, process, "bkg")
                elif process == "other":
                    if channel == "em":
                        other_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "W").Clone(), xlow=30.)
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "VVL"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ZL"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ggH125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "qqH125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "bbH125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ggHWW125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "qqHWW125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ZHWW125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "WHWW125"), xlow=30.))
                        plot.add_hist(other_hist, process, "bkg")
                    elif channel in ["et", "mt", "lt", "tt"]:
                        other_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "VVL").Clone(), xlow=30.)
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ZL"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ggH125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "qqH125"), xlow=30.))
                        other_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "bbH125"), xlow=30.))
                        plot.add_hist(other_hist, process, "bkg")
                else:
                    #print era, channel, category, process
                    plot.add_hist(
                        rebin_hist_for_logX(rootfile.get(era, channel, category, process), xlow=30.), process, "bkg")
                plot.setGraphStyle(
                    process, "hist", fillcolor=styles.color_dict[process])

            # get signal histograms
            if int(category) > 30:
                plot_idx_to_add_signal = [0,2] if args.linear else [1,2]
                for i in plot_idx_to_add_signal:
                    if args.model_independent:
                        ggH_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "ggh_t"), xlow=30.).Clone()
                        ggH_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ggh_i"), xlow=30.))
                        ggH_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "ggh_b"), xlow=30.))
                        plot.subplot(i).add_hist(ggH_hist, "ggH")
                        # bbH signal
                        bbH_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "bbh"), xlow=30.)
                        # if args.cross_section_bbh != args.cross_section_ggh:
                        #     print("Scaling bbH by {}".format(float(args.cross_section_bbh)/float(args.cross_section_ggh)))
                        #     bbH_hist.Scale(float(args.cross_section_bbh)/float(args.cross_section_ggh))
                        if args.cross_section_bbh != args.cross_section_ggh:
                            print("Scaling bbH by {}".format(float(args.cross_section_bbh)/0.0030584))
                            bbH_hist.Scale(float(args.cross_section_bbh)/0.0030584)
                        plot.subplot(i).add_hist(bbH_hist, "bbH")
                        # vector leptoquark signal
                        VLQ_hist = rebin_hist_for_logX(rootfile.get(era, channel, category, "VLQ_s"), xlow=30.).Clone()
                        VLQ_hist.Add(rebin_hist_for_logX(rootfile.get(era, channel, category, "VLQ_i"), xlow=30.))
                        plot.subplot(i).add_hist(
                            VLQ_hist, "VLQ")
                    else:
                        plot.subplot(i).add_hist(
                            rebin_hist_for_logX(rootfile.get(era, channel, category, "TotalSig"), xlow=30.), "mssm_sig")

            # get observed data and total background histograms
            plot.add_hist(
                rebin_hist_for_logX(rootfile.get(era, channel, category, "data_obs"), xlow=30.), "data_obs")
            plot.add_hist(
                rebin_hist_for_logX(rootfile.get(era, channel, category, "TotalBkg"), xlow=30.), "total_bkg")

            plot.subplot(0).setGraphStyle("data_obs", "e0")
            if args.model_independent:
                if int(category) > 30:
                    plot.subplot(0 if args.linear else 1).setGraphStyle(
                        "ggH", "hist", linecolor=styles.color_dict["ggH"], linewidth=2)
                    plot.subplot(0 if args.linear else 1).setGraphStyle(
                        "bbH", "hist", linecolor=styles.color_dict["bbH"], linewidth=2)
                    plot.subplot(0 if args.linear else 1).setGraphStyle(
                        "VLQ", "hist", linecolor=styles.color_dict["VLQ"], linewidth=2)
            else:
                plot.subplot(0 if args.linear else 1).setGraphStyle(
                    "mssm_sig", "hist", linecolor=styles.color_dict["bbH"], linewidth=2)
            plot.setGraphStyle(
                "total_bkg",
                "e2",
                markersize=0,
                fillcolor=styles.color_dict["unc"],
                linecolor=0)

            # assemble ratio
            if args.model_independent:
                if args.b_only_fit is not None:
                    rootfile_b_only = rootfile_parser.Rootfile_parser(args.b_only_fit, mode="CombineHarvesterMerged")
                    plot.subplot(2).add_hist(
                        rebin_hist_for_logX(rootfile_b_only.get(era, channel, category, "TotalBkg"), xlow=30.), "total_bkg_b_only")
                    plot.subplot(2).setGraphStyle(
                        "total_bkg_b_only", "hist", linecolor=ROOT.kBlue,
                        linewidth=2, linestyle=2)
                if int(category) > 30:
                    bkg_ggH = plot.subplot(2).get_hist("ggH")
                    bkg_bbH = plot.subplot(2).get_hist("bbH")
                    bkg_VLQ = plot.subplot(2).get_hist("VLQ")
                    bkg_ggH.Add(plot.subplot(2).get_hist("total_bkg"))
                    bkg_bbH.Add(plot.subplot(2).get_hist("total_bkg"))
                    bkg_VLQ.Add(plot.subplot(2).get_hist("total_bkg"))
                    plot.subplot(2).add_hist(bkg_ggH, "bkg_ggH")
                    plot.subplot(2).add_hist(bkg_bbH, "bkg_bbH")
                    plot.subplot(2).add_hist(bkg_VLQ, "bkg_VLQ")
                    plot.subplot(2).setGraphStyle(
                        "bkg_ggH",
                        "hist",
                        linecolor=styles.color_dict["ggH"],
                        linewidth=2)
                    plot.subplot(2).setGraphStyle(
                        "bkg_bbH",
                        "hist",
                        linecolor=styles.color_dict["bbH"],
                        linewidth=2)
                    plot.subplot(2).setGraphStyle(
                        "bkg_VLQ",
                        "hist",
                        linecolor=styles.color_dict["VLQ"],
                        linewidth=2)
                    to_normalize = ["total_bkg", "bkg_ggH", "bkg_bbH", "bkg_VLQ",
                        "data_obs"
                    ]
                    if args.b_only_fit is not None:
                        to_normalize.append("total_bkg_b_only")
                    plot.subplot(2).normalize(to_normalize, "total_bkg")
                else:
                    to_normalize = ["total_bkg", "data_obs" ]
                    if args.b_only_fit is not None:
                        to_normalize.append("total_bkg_b_only")
                    plot.subplot(2).normalize(to_normalize, "total_bkg")
            else:
                bkg_sig = plot.subplot(2).get_hist("mssm_sig")
                bkg_sig.Add(plot.subplot(2).get_hist("total_bkg"))
                plot.subplot(2).add_hist(bkg_sig, "bkg_mssm_sig")
                plot.subplot(2).add_hist(bkg_sig, "bkg_mssm_sig_top")
                plot.subplot(2).setGraphStyle(
                    "bkg_mssm_sig",
                    "hist",
                    linecolor=styles.color_dict["bbH"],
                    linewidth=3)
                plot.subplot(2).setGraphStyle(
                    "bkg_mssm_sig_top",
                    "hist",
                    linecolor=0)
                plot.subplot(2).normalize([
                    "total_bkg", "bkg_mssm_sig", "bkg_mssm_sig_top",
                    "data_obs"
                ], "total_bkg")

            # stack background processes
            plot.create_stack(bkg_processes, "stack")

            # Use binning from data histogram to get bin widths for the normalization.
            hist_for_rebinning = rootfile.get(era, channel, category, "data_obs")
            widths = []
            for i in range(hist_for_rebinning.GetNbinsX()):
                widths.append(hist_for_rebinning.GetBinWidth(i+1))
            # normalize stacks by bin-width
            if args.normalize_by_bin_width:
                plot.subplot(0).normalizeByBinWidth(widths=widths)
                plot.subplot(1).normalizeByBinWidth(widths=widths)

            if args.x_range is not None:
                for i in range(3):
                    plot.subplot(i).setXlims(*args.x_range)

            # set axes limits and labels
            if args.x_range is not None:
                range_hist = plot.subplot(0).get_hist("data_obs").Clone()
                range_hist.GetXaxis().SetRangeUser(*args.x_range)
                plot.subplot(0).setYlims(
                    split_dict[channel],
                    max(1.8 * range_hist.GetMaximum(),
                        split_dict[channel] * 2))
                    # max(1,
                    #     split_dict[channel] * 2))
            else:
                plot.subplot(0).setYlims(
                    split_dict[channel],
                    max(1.4 * plot.subplot(0).get_hist("data_obs").GetMaximum(),
                        split_dict[channel] * 2))

            if channel == "em":
                if int(category) == 36:
                    plot.subplot(2).setYlims(0.2, 2.5)
                else:
                    plot.subplot(2).setYlims(0.5, 2.0)
            elif channel == "tt":
                plot.subplot(2).setYlims(0.7, 2.5)
            else:
                if int(category) == 35:
                    plot.subplot(2).setYlims(0.5, 2.2)
                else:
                    plot.subplot(2).setYlims(0.7, 2.2)

            ROOT.TGaxis.SetMaxDigits(4)

            if args.linear != True:
                # plot.subplot(1).setYlims(1.e-4, split_dict[channel])
                plot.subplot(1).setYlims(min(1.e-3, plot.subplot(1).get_hist("data_obs").GetMinimum()/10., plot.subplot(1).get_hist("total_bkg").GetMinimum()/10.), split_dict[channel])
                plot.subplot(1).setLogY()
                if int(category) > 30 or int(category) == 2:
                    plot.subplot(1).setLogX()
                plot.subplot(1).setYlabel(
                    "")  # otherwise number labels are not drawn on axis
            if args.control_variable != None:
                if args.control_variable in styles.x_label_dict[args.channels[0]]:
                    x_label = styles.x_label_dict[args.channels[0]][
                        args.control_variable]
                else:
                    x_label = args.control_variable
                plot.subplot(2).setXlabel(x_label)
            else:
                if int(category) > 30 or int(category) == 2:
                    plot.subplot(2).setXlabel("m_{T}^{tot} (GeV)")
                else:
                    plot.subplot(2).setXlabel("SVFit m_{#tau#tau} (GeV)")
            if args.normalize_by_bin_width:
                if int(category) > 30 or int(category) == 2:
                    plot.subplot(0).setYlabel("dN/dm_{T}^{tot} (1/GeV)")
                    # plot.subplot(0).setYlabel("< Events / GeV >")
                else:
                    plot.subplot(0).setYlabel("dN/dm_{#tau#tau} (1/GeV)")
            else:
                plot.subplot(0).setYlabel("Events / {} GeV".format(plot.subplot(0).get_hist("total_bkg").GetBinWidth(1)))

            plot.subplot(2).setYlabel("Obs./Exp.")


            if (int(category) > 30 or int(category) == 2) and args.x_range is None:
                plot.subplot(0).setLogX()
                plot.subplot(2).setLogX()

            #plot.scaleXTitleSize(0.8)
            #plot.scaleXLabelSize(0.8)
            #plot.scaleYTitleSize(0.8)
            plot.scaleYLabelSize(0.8)
            #plot.scaleXLabelOffset(2.0)

            plot.scaleYTitleOffset(1.05)

            if channel == "em":
                if category == "32":
                    plot.subplot(0).setNYdivisions(3, 5)
                elif category == "33":
                    plot.subplot(0).setNYdivisions(4, 5)
            if int(category) < 35 and int(category) > 10:
                plot.subplot(1).setNYdivisions(3, 5)
            plot.subplot(2).setNYdivisions(3, 5)

            # draw subplots. Argument contains names of objects to be drawn in corresponding order.
            # procs_to_draw = ["stack", "total_bkg", "ggH", "ggH_top", "bbH", "bbH_top", "data_obs"] if args.linear else ["stack", "total_bkg", "data_obs"]
            if category == "2":
                procs_to_draw = ["stack", "total_bkg", "data_obs"] if args.linear else ["stack", "total_bkg", "data_obs"]
            else:
                if args.model_independent:
                    procs_to_draw = ["stack", "total_bkg", "ggH", "bbH", "data_obs"] if args.linear else ["stack", "total_bkg", "data_obs"]
                else:
                    procs_to_draw = ["stack", "total_bkg", "mssm_sig", "data_obs"] if args.linear else ["stack", "total_bkg", "data_obs"]
                if args.blinded:
                    procs_to_draw.remove("data_obs")
            plot.subplot(0).Draw(procs_to_draw)
            if not args.linear:
                # plot.subplot(1).Draw([
                #     "stack", "total_bkg", "ggH", "bbH",
                #     "ggH_top", "bbH_top",
                #     "data_obs"
                # ])
                if category == "2":
                    plot.subplot(1).Draw([
                        "stack", "total_bkg",
                        "data_obs"
                    ])
                else:
                    if args.model_independent:
                        procs_to_draw = ["stack", "total_bkg", "ggH", "bbH", "VLQ", "data_obs"]
                    else:
                        procs_to_draw = ["stack", "total_bkg", "mssm_sig", "data_obs"]
                    if args.blinded:
                        procs_to_draw.remove("data_obs")
                    plot.subplot(1).Draw(procs_to_draw)

                plot.subplot(0).remove_lower_x_ticks()
                plot.subplot(1)._pad.SetTickx(0)
                plot.add_line(0, 30, split_dict[channel], 6000, split_dict[channel], linestyle=1, color=ROOT.kGray+2)
                for line in plot._lines:
                    line.Draw()
                label2 = ROOT.TLatex()
                label2.SetNDC()
                label2.SetTextAngle(270)
                label2.SetTextColor(ROOT.kGray+2)
                label2.SetTextSize(0.030)
                label2.DrawLatex(0.97, 0.54, "log scale")
                label2.DrawLatex(0.97, 0.75, "linear scale")
                # Redraw upper x axis to have good looking plot
                plot.subplot(0)._pad.cd()
                y_max = 1.4 * plot.subplot(0).get_hist("data_obs").GetMaximum()
                axis = ROOT.TGaxis(30., y_max, 5000., y_max, 30., 5000., 510, "GBSU-")
                axis.SetTickLength(0.02)
                axis.SetTickSize(0.02)
                axis.Draw()

            # Draw upper panel after lower panel of split to ensure full
            # display of the data points
            if category == "2":
                plot.subplot(2).Draw([
                    "total_bkg",
                    "data_obs"
                ])
            else:
                if args.model_independent:
                    if args.b_only_fit is not None:
                        if int(category) > 30:
                            procs_to_draw = ["total_bkg_b_only", "total_bkg", "bkg_ggH", "bkg_bbH", "bkg_VLQ", "data_obs"]
                        else:
                            procs_to_draw = ["total_bkg_b_only", "total_bkg", "data_obs"]
                    else:
                        if int(category) > 30:
                            procs_to_draw = ["total_bkg", "bkg_ggH", "bkg_bbH", "bkg_VLQ", "data_obs"]
                        else:
                            procs_to_draw = ["total_bkg", "data_obs"]
                else:
                    procs_to_draw = ["total_bkg", "bkg_mssm_sig", "data_obs"]
                if args.blinded:
                    procs_to_draw.remove("data_obs")
                plot.subplot(2).Draw(procs_to_draw)

            # create legends
            if int(category) < 30 and int(category) > 2:
                plot.add_legend(width=0.38, height=0.30)
            else:
                # plot.add_legend(width=0.60, height=0.20, pos=1)
                plot.add_legend(width=0.40, height=0.33, pos=3)
            # plot.add_legend(width=0.6, height=0.15)
            for process in legend_bkg_processes:
                if process == "HSM" and channel == "em":
                    plot.legend(0).add_entry(
                        0, process, "H(125 GeV)", "f")
                else:
                    plot.legend(0).add_entry(
                        0, process, styles.legend_label_dict[process.replace("TTL", "TT").replace("VVL", "VV")], 'f')
            if not args.blinded:
                plot.legend(0).add_entry(0, "data_obs", "Observed", 'PE')
            plot.legend(0).add_entry(0, "total_bkg", "Bkg. unc.", 'f')
            if args.control_region and category == "2":
                pass
            else:
                if args.model_independent:
                    if int(category) > 30:
                        mass_str = "%s GeV" % args.mass if float(args.mass) < 1000 else "{:.1f} TeV".format(float(args.mass) / 1000)
                        ggH_xs_str = "%s pb" % args.cross_section_ggh if float(args.cross_section_ggh) > 1e-2 else "{} fb".format(float(args.cross_section_ggh) * 1000)
                        bbH_xs_str = "%s pb" % args.cross_section_bbh if float(args.cross_section_bbh) > 1e-2 else "{} fb".format(float(args.cross_section_bbh) * 1000)
                        plot.legend(0).add_entry(0 if args.linear else 1, "ggH", "#splitline{gg#phi @ %s}{(m_{#phi} = %s)}" % (ggH_xs_str, mass_str), 'l')
                        plot.legend(0).add_entry(0 if args.linear else 1, "bbH", "#splitline{bb#phi @ %s}{(m_{#phi} = %s)}" % (bbH_xs_str, mass_str), 'l')
                        plot.legend(0).add_entry(0 if args.linear else 1, "", "", '')
                        plot.legend(0).add_entry(0 if args.linear else 1, "VLQ", "#splitline{VLQ, g_{U} = %s}{(m_{U} = %s TeV)}" % (1.2, 1), 'l')
                else:
                    plot.legend(0).add_entry(0 if args.linear else 1, "mssm_sig", "#splitline{H #rightarrow #tau#tau}{#splitline{(m_{A}= %s GeV,}{ tan #beta = %s)}}" %(args.mass, args.tanbeta), 'l')
            plot.legend(0).setNColumns(2)
            # plot.legend(0).scaleTextSize(1.05)
            plot.legend(0).Draw()


            # FIXME: Legend for ratio plot temporarily disabled.
            if channel != "em":
                plot.add_legend(
                    reference_subplot=2, pos=1, width=0.45, height=0.10)
            else:
                plot.add_legend(
                    reference_subplot=2, pos=1, width=0.75, height=0.05)
                if not args.blinded:
                    plot.legend(1).add_entry(0, "data_obs", "Observed", 'PE')
                plot.legend(1).add_entry(0, "total_bkg", "Bkg. unc.", 'f')
                if args.b_only_fit:
                    plot.legend(1).add_entry(2, "total_bkg_b_only", "Bkg. only fit", "l")
            if args.control_region and category == "2":
                pass
            else:
                if args.model_independent:
                    if int(category) > 30:
                        plot.legend(1).add_entry(0 if args.linear else 1, "ggH",
                                             "gg#phi", 'l')
                        plot.legend(1).add_entry(0 if args.linear else 1, "bbH",
                                             "bb#phi", 'l')
                        plot.legend(1).add_entry(0 if args.linear else 1, "VLQ",
                                             "VLQ", 'l')
                else:
                    plot.legend(1).add_entry(0 if args.linear else 1, "mssm_sig",
                                                 "H+bkg.", 'l')
            if channel != "em":
                if not args.blinded:
                    plot.legend(1).add_entry(0, "data_obs", "Observed", 'PE')
                plot.legend(1).add_entry(0, "total_bkg", "Bkg. unc.", 'f')
                if args.b_only_fit is not None:
                    plot.legend(1).add_entry(2, "total_bkg_b_only", "Bkg. only fit", "l")
                plot.legend(1).setNColumns(3)
            else:
                plot.legend(1).setNColumns(5)
                if args.b_only_fit is not None:
                    plot.legend(1).setNColumns(6)
            # plot.legend(1).scaleTextSize(1.05)
            plot.legend(1).Draw()

            # draw additional labels
            plot.DrawCMS(cms_sub="")
            if "2016" in args.era:
                plot.DrawLumi("36.3 fb^{-1} (2016, 13 TeV)")
            elif "2017" in args.era:
                plot.DrawLumi("41.5 fb^{-1} (2017, 13 TeV)")
            elif "2018" in args.era:
                plot.DrawLumi("59.7 fb^{-1} (2018, 13 TeV)")
            elif "combined" in args.era:
                plot.DrawLumi("138 fb^{-1} (13 TeV)")
            else:
                logger.critical("Era {} is not implemented.".format(args.era))
                raise Exception

            posChannelCategoryLabelLeft = None
            plot.DrawChannelCategoryLabel(
                "#font[42]{%s, %s}" % (channel_dict[channel], category_dict[channel][category]),
                begin_left=posChannelCategoryLabelLeft)

            # save plot
            if not os.path.exists(args.output_dir):
                os.makedirs(args.output_dir)
            postfix = "prefit" if "prefit" in args.input else "postfit" if "postfit" in args.input else "undefined"
            plot.save("%s/%s_%s_%s_%s.%s" % (args.output_dir, args.era, channel, args.control_variable if args.control_variable is not None else category,
                                                postfix, "png"
                                                if args.png else "pdf"))
            # Undo switch of split value
            if int(category) < 35 and int(category) > 10:
                split_dict[channel] -= 100
            plots.append(
                plot
            )  # work around to have clean up seg faults only at the end of the script


if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("{}_plot_shapes.log".format(args.era), logging.INFO)
    main(args)
