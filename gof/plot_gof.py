#!/usr/bin/env python

import argparse
import os
import json
import yaml
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['font.size'] = 16
import matplotlib.pyplot as plt
from matplotlib import cm

import logging
logger = logging.getLogger("plot_gof")
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

category_dict = {
    "et": {
        "32": r"No b-tag, Tight $m_{T}$",
        "33": r"No b-tag, Loose $m_{T}$",
        "35": r"b-tag, Tight $m_{T}$",
        "36": r"b-tag, Loose $m_{T}$",
    },
    "mt": {
        "32": r"No b-tag, Tight $m_{T}$",
        "33": r"No b-tag, Loose $m_{T}$",
        "35": r"b-tag, Tight $m_{T}$",
        "36": r"b-tag, Loose $m_{T}$",
    },
    "em": {
        "2": "ttbar control",
        "32": r"No b-tag, high $d_{\zeta}$",
        "33": r"No b-tag, medium $d_{\zeta}$",
        "34": r"No b-tag, low $d_{\zeta}$",
        "35": r"b-tag, high $d_{\zeta}$",
        "36": r"b-tag, medium $d_{\zeta}$",
        "37": r"b-tag, low $d_{\zeta}$",
    },
    "tt": {
        "32": "No b-tag",
        "35": "b-tag",
    },
}

channel_dict = {
    "et": r"$e\tau_{h}$",
    "mt": r"$\mu\tau_{h}$",
    "tt": r"$\tau_{h}\tau_{h}$",
    "em": r"$e\mu$",
}

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot goodness of fit results")
    parser.add_argument(
        "path", help="Path to directory with goodness of fit results")
    parser.add_argument("channel", type=str, help="Select channel to be plotted")
    parser.add_argument("era", type=str, help="Select era to be plotted")
    parser.add_argument("mode", choices=["cats-per-era", "cats-per-channel", "per-channel"],
                        help="Plotting mode.")
    parser.add_argument("-c", "--comparison-path",
                        default=None,
                        help="Path to directory with GoF results for comparison.")
    return parser.parse_args()


def make_cmap(colors, position):
    if len(position) != len(colors):
        raise Exception("position length must be the same as colors")
    elif position[0] != 0 or position[-1] != 1:
        raise Exception("position must start with 0 and end with 1")

    cdict = {'red': [], 'green': [], 'blue': []}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 256)
    return cmap


def plot_1d(categories, results, filename, channel, comparison=None):
    plt.figure(figsize=(len(categories) * 0.5, 5.0))
    y = results
    x = range(len(y))
    plt.plot(x, y, '+', mew=4, ms=16)
    if comparison is not None:
        plt.plot(x, comparison, '+', mew=4, ms=16)
    plt.ylim((-0.05, 1.05))
    plt.xlim((-0.5, len(x) - 0.5))
    plt.xticks(x, categories, rotation='vertical')
    plt.axhline(y=0.05, linewidth=3, color='r')
    for i, res in enumerate(y):
        if res < 0.05:
            plt.text(i, 0.9, "{:.3f}".format(res), rotation='vertical',
                     horizontalalignment="center", verticalalignment="center")
    plt.ylabel('Saturated goodness of fit p-value', labelpad=20)
    ax = plt.gca()
    ax.xaxis.grid()
    plt.savefig(filename+".png" if comparison is None else filename+"-comparison.png", bbox_inches="tight")
    plt.savefig(filename+".pdf" if comparison is None else filename+"-comparison.pdf", bbox_inches="tight")


def search_results_1d(path, mode, era="combined", channel="cmb", forComp=False):
    results = []
    missing = []
    if mode == "cats-per-era":
        for ch_ in ["et", "mt", "tt", "em"]:
            for category in sorted(category_dict[ch_].keys()):
                p_val, present = get_gof_result(path, era, ch_, category, forComp)
                if present:
                    results.append(p_val)
                else:
                    missing.append("-".join([era, ch_, category]))
    elif mode == "cats-per-channel":
        for era_ in ["2016", "2017", "2018"]:
            for category in sorted(category_dict[channel].keys()):
                p_val, present = get_gof_result(path, era_, channel, category)
                if present:
                    results.append(p_val)
                else:
                    missing.append("-".join([era_, channel, category]))
    elif mode == "per-channel":
        for era_ in ["2016", "2017", "2018"]:
            for ch_ in ["et", "mt", "tt", "em"]:
                p_val, present = get_gof_result(path, era_, ch_)
                if present:
                    results.append(p_val)
                else:
                    missing.append("-".join([era_, ch_]))
    return missing, results


def get_gof_result(path, era, channel, category=None, forComp=False):
    if category is None:
        filename = os.path.join(path, era, channel, "gof",
                                "{}-{}".format(era, channel), "gof.json" if not forComp else "gof.json")
    else:
        filename = os.path.join(path, era, channel, "gof",
                                "{}-{}-{}".format(era, channel, category), "gof.json" if not forComp else "gof.json")
    print("Looking for file with name: {}".format(filename))
    if not os.path.exists(filename):
        return -1.0, False
    logger.debug(
        "Found goodness of fit result for category %s in channel %s in era %s.",
        category, channel, era)

    p_value = json.load(open(filename))
    return p_value["160.0"]["p"], True


def main(args):
    if not os.path.exists(args.path):
        logger.fatal("Path %s does not exist.")
        raise Exception

    # Plot 1D gof results
    missing_1d, results_1d = search_results_1d(args.path, args.mode, args.era, args.channel)
    logger.debug("Missing categories for 1D plot in channel %s:", args.channel)
    for category in missing_1d:
        print("{}".format(category))
    if args.comparison_path is not None:
        missing_1d_comp, results_1d_comp = search_results_1d(args.comparison_path, args.mode, args.era, args.channel, forComp=True)
        logger.debug("Missing categories for 1D comparison plot in channel %s:", args.channel)
        for category in missing_1d_comp:
            print("{}".format(category))

    avail_era = ["2016", "2017", "2018"]
    avail_ch = ["et", "mt", "tt", "em"]
    if args.mode == "cats-per-era":
        categories = [", ".join([channel_dict[ch_], category_dict[ch_][x]]) for ch_ in avail_ch for x in sorted(category_dict[ch_].keys())]
        outname = "{}/{}_gof_{}".format(args.path, args.era, args.mode)
    elif args.mode == "cats-per-channel":
        categories = [", ".join([era_, category_dict[args.channel][x]]) for era_ in avail_era for x in sorted(category_dict[args.channel].keys())]
        outname = "{}/{}_gof_{}".format(args.path, args.channel, args.mode)
    elif args.mode == "per-channel":
        categories = [", ".join([era_, channel_dict[ch_]]) for era_ in avail_era for ch_ in avail_ch]
        outname = "{}/gof_{}".format(args.path, args.mode)

    if args.comparison_path is None:
        plot_1d(categories, results_1d, outname, args.channel)
    else:
        plot_1d(categories, results_1d, outname, args.channel, comparison=results_1d_comp)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
