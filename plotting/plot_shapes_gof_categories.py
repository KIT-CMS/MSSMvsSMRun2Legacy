#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Dumbledraw.dumbledraw as dd
import Dumbledraw.rootfile_parser as rootfile_parser
import Dumbledraw.styles as styles
import ROOT

import argparse
import copy
import yaml
import distutils.util
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
    parser.add_argument("-o", "--outputfolder", type=str, required=True, help="...yourself")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="ROOT file with shapes of processes")
    parser.add_argument(
        "--gof-variable",
        type=str,
        default=None,
        help="Enable plotting goodness of fit shapes for given variable")
    parser.add_argument(
        "--png", action="store_true", help="Save plots in png format")
    parser.add_argument(
        "--categories",
        type=str,
        required=True,
        help="Select categorization.")
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
        "--background-only",
        type=lambda x:bool(distutils.util.strtobool(x)),
        default=False,
        help="Plot only the background categories")
    parser.add_argument(
        "--chi2test",
        action="store_true",
        help="Print chi2/ndf result in upper-right of subplot")

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


def check_for_zero_bins(hist):
    for ibin in range(1, hist.GetNbinsX()+1):
        if hist.GetBinContent(ibin) < 1.e-9:
            logger.debug("Found bin with content %d in histogram %s. "
                         "Setting bin content to zero.",
                         hist.GetBinContent(ibin),
                         hist.GetName())
            hist.SetBinContent(ibin, 0)
    return hist


def get_histogram(rootfile, era, channel, category, process):
    if era == "combined":
        hist = get_histogram(rootfile, "2016", channel, category, process).Clone()
        hist.Add(get_histogram(rootfile, "2017", channel, category, process))
        hist.Add(get_histogram(rootfile, "2018", channel, category, process))
    else:
        if channel == "tt" and process == "jetFakes":
            hist = rootfile.get(era, channel, category, process).Clone()
            hist.Sumw2()
            hist.Add(rootfile.get(era, channel, category, "wFakes"))
        elif channel == "em" and process == "EWK":
            hist = rootfile.get(era, channel, category, "W").Clone()
            hist.Sumw2()
            hist.Add(rootfile.get(era, channel, category, "VVL"))
        else:
            hist = rootfile.get(era, channel, category, process)
            hist.Sumw2()
    return hist


def main(args):
    #### plot signals
    if args.background_only:
        stxs_stage1p1_cats = []
    else:
        stxs_stage1p1_cats = [str(100+i) for i in range(5)] + [str(200+i) for i in range(4)]
    print(args)
    if args.gof_variable != None:
        channel_categories = {
            "et": ["300"],
            "mt": ["301", "302"],
            "tt": ["300"],
            "em": ["301", "302"],
        }
    else:
        channel_categories = {
            #"et": ["ztt", "zll", "w", "tt", "ss", "misc"],
            "et": ["12", "15", "11", "13", "14", "16"],
            #"mt": ["ztt", "zll", "w", "tt", "ss", "misc"],
            "mt": ["12", "15", "11", "13", "14", "16"],
            #"tt": ["ztt", "noniso", "misc"]
            "tt": ["12", "17", "16"],
            #"em": ["ztt", "tt", "ss", "misc", "db"]
            "em": ["12", "13", "14", "16", "19"]
        }

    channel_dict = {
        "ee": "ee",
        "em": "e#mu",
        "et": "e#tau_{h}",
        "mm": "#mu#mu",
        "mt": "#mu#tau_{h}",
        "tt": "#tau_{h}#tau_{h}"
    }
    if args.gof_variable != None:
        category_dict = {
                "300": "inclusive",
                "301": "No b-tag",
                "302": "b-tag",
        }
    else:
        category_dict = {
            "1": "ggh",
            "100": "ggh 0-jet",
            "101": "ggh 1-jet p_{T}^{H} [0,120]",
            "102": "ggh 1-jet p_{T}^{H} [120,200]",
            "103": "ggh #geq 2-jet",
            "104": "ggh p_{T}^{H}>200",
            "2": "qqh",
            "200": "qqh 2J low mjj",
            "201": "qqh p_{T}^{H}>200",
            "202": "qqh vbftopo mjj>700",
            "203": "qqh vbftopo mjj [350,700]",
            "12": "ztt",
            "15": "zll",
            "11": "wjets",
            "13": "tt",
            "14": "qcd",
            "16": "misc",
            "17": "qcd",
            "19": "diboson",
            "20": "Genuine #tau",
            "21": "Jet #rightarrow #tau_{h}"
        }
    if args.linear == True:
        split_value = 0
    else:
        if args.normalize_by_bin_width:
            split_value = 10001
        else:
            split_value = 101

    split_dict = {c: split_value for c in ["et", "mt", "tt", "em"]}

    bkg_processes = [
        "ggH125", "qqH125", "VVL", "TTL", "ZL", "jetFakes", "EMB"
    ]
    if not args.fake_factor and args.embedding:
        bkg_processes = [
            "QCD", "VVJ", "W", "TTJ", "ZJ", "ZL", "EMB"
        ]
    if not args.embedding and args.fake_factor:
        bkg_processes = [
            "VVT", "VVL", "TTT", "TTL", "ZL", "jetFakes", "ZTT"
        ]
    if not args.embedding and not args.fake_factor:
        bkg_processes = [
            "QCD", "VVT", "VVL", "VVJ", "W", "TTT", "TTL", "TTJ", "ZJ", "ZL", "ZTT"
        ]
    all_bkg_processes = [b for b in bkg_processes]
    legend_bkg_processes = copy.deepcopy(bkg_processes)
    legend_bkg_processes.reverse()

    if "2016" in args.era:
        era = "2016"
    elif "2017" in args.era:
        era = "2017"
    elif "2018" in args.era:
        era = "2018"
    elif "combined" in args.era:
        era = "combined"
    else:
        logger.critical("Era {} is not implemented.".format(args.era))
        raise Exception
    print(channel_categories)
    plots = []
    for channel in args.channels:
        for category in channel_categories[channel]:
            print "Plot for category: ",category
            rootfile = rootfile_parser.Rootfile_parser(args.input)
            if channel == "em" and args.embedding:
                bkg_processes = ["ggH125", "qqH125", "TTL", "QCD", "EWK", "ZL", "EMB"]
            elif channel == "em" and not args.embedding:
                bkg_processes = ["VVL", "W", "TTL", "ZL", "QCD", "ZTT"]
            else:
                bkg_processes = [b for b in all_bkg_processes]
            legend_bkg_processes = copy.deepcopy(bkg_processes)
            legend_bkg_processes.reverse()
            # create plot
            width = 600
            if args.linear == True:
                plot = dd.Plot(
                    [0.25, [0.25, 0.23]], "ModTDR", r=0.04, l=0.14, width=width)
            else:
                plot = dd.Plot(
                    [0.5, [0.3, 0.28]], "ModTDR", r=0.04, l=0.14, width=width)

            # get background histograms
            for process in bkg_processes:
                try:
                    plot.add_hist(get_histogram(rootfile, era, channel, category, process),
                                  process, "bkg")
                    plot.setGraphStyle(
                        process, "hist", fillcolor=styles.color_dict[process])
                except:
                    pass


            # get observed data and total background histograms
            # NOTE: With CMSSW_8_1_0 the TotalBkg definition has changed.
            print("plot.add_hist(rootfile.get("+era+", "+channel+", "+category+', "data_obs")')
            plot.add_hist(
                get_histogram(rootfile, era, channel, category, "data_obs"), "data_obs")
            total_bkg = check_for_zero_bins(get_histogram(rootfile, era, channel, category, "TotalBkg"))
            #ggHHist = rootfile.get(era, channel, category, "ggH")
            #qqHHist = rootfile.get(era, channel, category, "qqH")
            #total_bkg.Add(ggHHist, -1)
            #if qqHHist:
            #     total_bkg.Add(qqHHist, -1)
            plot.add_hist(total_bkg, "total_bkg")

            plot.subplot(0).setGraphStyle("data_obs", "e0")
            plot.setGraphStyle(
                "total_bkg",
                "e2",
                markersize=0,
                fillcolor=styles.color_dict["unc"],
                linecolor=0)

            plot.subplot(2).normalize([
                "total_bkg",
                "data_obs"
            ], "total_bkg")

            # stack background processes
            plot.create_stack(bkg_processes, "stack")

            # normalize stacks by bin-width
            if args.normalize_by_bin_width:
                plot.subplot(0).normalizeByBinWidth()
                plot.subplot(1).normalizeByBinWidth()

            # set axes limits and labels
            if channel == "em":
                # plot.setXlims(-200.,150.)
                if category in ["300", "301"]:
                    plot.setXlims(-70., 100.)
                    plot.subplot(0).setYlims(
                        split_dict[channel],
                        max(1.40 * plot.subplot(0).get_hist("data_obs").GetMaximum(),
                            split_dict[channel] * 2))
                elif category == "302":
                    plot.setXlims(-140., 140)
                    plot.subplot(0).setYlims(
                        split_dict[channel],
                        max(1.40 * plot.subplot(0).get_hist("data_obs").GetMaximum(),
                            split_dict[channel] * 2))
            elif channel in ["et", "mt"]:
                if category in ["300", "301"]:
                    plot.subplot(0).setYlims(
                        split_dict[channel],
                        max(1.25 * plot.subplot(0).get_hist("data_obs").GetMaximum(),
                            split_dict[channel] * 2))
                elif category in ["302"]:
                    plot.subplot(0).setYlims(
                        split_dict[channel],
                        max(2.00 * plot.subplot(0).get_hist("data_obs").GetMaximum(),
                            split_dict[channel] * 2))

            if args.gof_variable is None:
                plot.subplot(2).setYlims(0.45, 2.05)
            else:
                plot.subplot(2).setYlims(0.85, 1.15)

            if not args.linear:
                plot.subplot(1).setYlims(0.1, split_dict[channel])
                plot.subplot(1).setLogY()
                plot.subplot(1).setYlabel(
                    "")  # otherwise number labels are not drawn on axis
            if args.gof_variable != None:
                gof_log_vars = []
                if args.gof_variable in gof_log_vars:
                    plot.subplot(0).setLogX()
                    plot.subplot(1).setLogX()
                    plot.subplot(2).setLogX()
                elif args.gof_variable in ["njets_red", "nbtag_red"]:
                    bins = plot.subplot(2).get_hist("data_obs").GetNbinsX()
                    bin_labels = map(str,range(bins-1))
                    bin_labels.append("#geq%i" % (bins-1) )
                    plot.subplot(2).setNXdivisions(bins, 0, 2)
                    plot.subplot(2).changeXLabels(bin_labels)
                if args.gof_variable in styles.x_label_dict[args.channels[0]]:
                    x_label = styles.x_label_dict[args.channels[0]][
                        args.gof_variable]
                else:
                    x_label = args.gof_variable
                plot.subplot(2).setXlabel(x_label)
            elif args.gof_variable != None and args.linear:
                if args.gof_variable in styles.x_label_dict[args.channels[0]]:
                    x_label = styles.x_label_dict[args.channels[0]][
                        args.gof_variable]
                else:
                    x_label = args.gof_variable
                plot.subplot(2).setXlabel(x_label)
            else:
                plot.subplot(2).setXlabel("NN output")
            if args.normalize_by_bin_width:
                plot.subplot(0).setYlabel("dN/d(%s)"%args.gof_variable)
            else:
                plot.subplot(0).setYlabel("Events / 10 GeV")

            plot.subplot(2).setYlabel("Data/Exp")
            plot.subplot(2).setGrid()

            #plot.scaleXTitleSize(0.8)
            #plot.scaleXLabelSize(0.8)
            #plot.scaleYTitleSize(0.8)
            plot.scaleYLabelSize(0.8)
            #plot.scaleXLabelOffset(2.0)
            plot.scaleYTitleOffset(1.1)

            plot.subplot(2).setNYdivisions(3, 0)

            #if not channel == "tt" and category in ["11", "12", "13", "14", "15", "16"]:
            #    plot.subplot(2).changeXLabels(["0.2", "0.4", "0.6", "0.8", "1.0"])

            # draw subplots. Argument contains names of objects to be drawn in corresponding order.
            procs_to_draw = ["stack", "total_bkg", "data_obs"] if args.linear else ["stack", "total_bkg", "data_obs"]
            plot.subplot(0).Draw(procs_to_draw)
            if args.linear != True:
                plot.subplot(1).Draw([
                    "stack", "total_bkg",
                    "data_obs"
                ])
            plot.subplot(2).Draw([
                "total_bkg",
                "data_obs"
            ])

            # create legends
            if channel == "mt" and category in ["302"]:
                plot.add_legend(width=0.5, height=0.20)
                plot.legend(0).setNColumns(2)
            elif channel == "em" and category in ["302"]:
                plot.add_legend(width=0.25, height=0.35)
                plot.legend(0).setNColumns(1)
            else:
                plot.add_legend(width=0.3, height=0.35)
                plot.legend(0).setNColumns(1)
            plot.legend(0).add_entry(0, "data_obs", "Data", 'PE')
            for process in legend_bkg_processes:
                try:
                    plot.legend(0).add_entry(
                        0, process, styles.legend_label_dict[process.replace("TTL", "TT").replace("VVL", "VV")], 'f')
                except:
                    pass
            plot.legend(0).add_entry(0, "total_bkg", "Bkg. unc.", 'f')

            plot.legend(0).Draw()

            if channel == "em":
                if category in ["300", "301"]:
                    _ymax_low = max(1.05 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    _ymax_high = max(0.6 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    plot.add_line(xmin=-35, xmax=-35, ymin=0, ymax=_ymax_low, linestyle=7, color=ROOT.kBlack)
                    plot.add_line(xmin=-10, xmax=-10, ymin=0, ymax=_ymax_low, linestyle=7, color=ROOT.kBlack)
                    plot.add_line(xmin= 30, xmax= 30, ymin=0, ymax=_ymax_low, linestyle=7, color=ROOT.kBlack)
                    plot.DrawText(0.335, 0.80, "low-p_{#zeta}", textsize=0.025)
                    plot.DrawText(0.465, 0.80, "medium-p_{#zeta}", textsize=0.025)
                    plot.DrawText(0.75, 0.5, "high-p_{#zeta}", textsize=0.025)
                elif category in ["302"]:
                    _ymax_low = max(1.05 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    _ymax_high = max(0.6 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    plot.add_line(xmin=-35, xmax=-35, ymin=0, ymax=_ymax_low, linestyle=7, color=ROOT.kBlack)
                    plot.add_line(xmin=-10, xmax=-10, ymin=0, ymax=_ymax_low, linestyle=7, color=ROOT.kBlack)
                    plot.add_line(xmin= 30, xmax= 30, ymin=0, ymax=_ymax_low, linestyle=7, color=ROOT.kBlack)
                    plot.DrawText(0.450, 0.8, "low-p_{#zeta}", textsize=0.025)
                    plot.DrawText(0.525, 0.8, "medium-p_{#zeta}", textsize=0.025)
                    plot.DrawText(0.75, 0.45, "high-p_{#zeta}", textsize=0.025)
                    plot.DrawText(0.305, 0.8, "t#bar{t} CR", textsize=0.025)
            elif channel in ["et", "mt"]:
                if category in ["300", "301"]:
                    _ymax = max(0.8 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    plot.add_line(xmin=40, xmax=40, ymin=0, ymax=_ymax, linestyle=7, color=ROOT.kBlack)
                    _ymax = max(0.5 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    plot.add_line(xmin=70, xmax=70, ymin=0, ymax=_ymax, linestyle=7, color=ROOT.kBlack)
                    plot.DrawText(0.24, 0.68, "Tight-m_{T}", textsize=0.030)
                    plot.DrawText(0.40, 0.55, "Loose-m_{T}", textsize=0.030)
                elif category in ["302"]:
                    _ymax = max(1.3 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    plot.add_line(xmin=40, xmax=40, ymin=0, ymax=_ymax, linestyle=7, color=ROOT.kBlack)
                    _ymax = max(1.3 * plot.subplot(0).get_hist("data_obs").GetMaximum(), split_dict[channel] * 2)
                    plot.add_line(xmin=70, xmax=70, ymin=0, ymax=_ymax, linestyle=7, color=ROOT.kBlack)
                    plot.DrawText(0.18, 0.67, "Tight-m_{T}", textsize=0.030)
                    plot.DrawText(0.40, 0.67, "Loose-m_{T}", textsize=0.030)

            for line in plot._lines:
                line.Draw()

            if args.chi2test:
                import ROOT as r
                f = r.TFile(args.input, "read")
                background = f.Get("htt_{}_{}_Run{}_{}/TotalBkg".format(
                    channel, category, args.era, "prefit"
                    if "prefit" in args.input else "postfit"))
                data = f.Get("htt_{}_{}_Run{}_{}/data_obs".format(
                    channel, category, args.era, "prefit"
                    if "prefit" in args.input else "postfit"))
                chi2 = data.Chi2Test(background, "UW CHI2/NDF")
                plot.DrawText(0.7, 0.3,
                              "\chi^{2}/ndf = " + str(round(chi2, 3)))

            plot.add_legend(
                reference_subplot=2, pos=1, width=0.5, height=0.03)
            plot.legend(1).add_entry(0, "data_obs", "Data", 'PE')
            plot.legend(1).add_entry(0, "total_bkg", "Bkg. unc.", 'f')
            plot.legend(1).setNColumns(4)
            plot.legend(1).Draw()

            # draw additional labels
            plot.DrawCMS()
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

            plot.DrawChannelCategoryLabel(
                "%s, %s" % (channel_dict[channel], category_dict[category]),
                begin_left=None)

            # save plot
            postfix = "prefit" if "prefit" in args.input else "postfit" if "postfit" in args.input else "undefined"
            plot.save("%s/%s_%s_%s_%s.%s" % (args.outputfolder, args.era, channel, category, postfix, "png" if args.png else "pdf"))
            plots.append(
                plot
            )  # work around to have clean up seg faults only at the end of the script


if __name__ == "__main__":
    args = parse_arguments()
    setup_logging("{}_plot_shapes.log".format(args.era), logging.INFO)
    main(args)
