from __future__ import print_function
import ROOT as r
import argparse
from prettytable import PrettyTable
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(
    description="Check Integrals of MSSM singals in SM categories")

parser.add_argument("--inputfile",
                    type=str,
                    required=True,
                    help="path to shapes")

parser.add_argument("--mssmfile",
                    type=str,
                    default="",
                    help="path to MSSM nominal signal shapes")

parser.add_argument('--channel',
                    type=str,
                    required=True,
                    help="Channel to analyze")
parser.add_argument('--era', type=str, required=True, help="Era to analyze")
parser.add_argument(
    '--threshold',
    type=float,
    default=5.0,
    help="threshold at which a process is considered significant")
parser.add_argument(
    '--mask',
    action='store_true',
    help="if set, all values below the threshold are masked from the plot")

parser.add_argument(
    '--entries',
    action='store_true',
    help="if set, the number of entries, instead of the integral is checked")

parser.add_argument(
    '--add-fraction',
    action='store_true',
    help="if set, the fraction of events ending up in the SM categories is shown")


args = parser.parse_args()
r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)
sns.set_style("ticks")
plt.rcParams.update({'font.size': 12})
plt.rc('ytick', labelsize=12)

category_dict = {
    "et": ["xxh", "tt", "zll", "misc", "emb", "ff"],
    "mt": ["xxh", "tt", "zll", "misc", "emb", "ff"],
    "em": ["xxh", "emb", "tt", "db", "misc", "ss"],
    "tt": ["xxh", "misc", "emb", "ff"],
}


def get_integrals(channel, f, mssmfile, categories, add_fraction, entries):
    output_data = {}
    output_data["mass"] = {}
    rfile = r.TFile.Open(f, "read")
    print(categories)
    print("Looking into", f)
    processes = sorted([
        p.GetName() for p in rfile.Get(categories[0]).GetListOfKeys()
        if "Up" not in p.GetName() and "Down" not in p.GetName()
        and "data_obs" not in p.GetName() and "2H" not in p.GetName()
    ])
    processes = [
        process for process in processes
        if any(signalprocess in process for signalprocess in [
            "ggh_t", "ggh_b", "ggh_i", "ggH_t", "ggH_b", "ggH_i", "ggA_t",
            "ggA_b", "ggA_i", "bbH"
        ])
    ]
    for category in categories:
        output_data[category] = {}
        print("\tCategory:", category)
        d = rfile.Get(category)
        for process in processes:
            mass = int(process.split("_")[-1])
            if process not in output_data["mass"]:
                output_data["mass"][process] = mass
            if args.entries:
                output_data[category][process] = d.Get(process).GetEntries()
            else:
                output_data[category][process] = d.Get(process).Integral()
    rfile.Close()

    if add_fraction:
        rmssmfile = r.TFile.Open(mssmfile, "read")
        fractions = {}
        for category in categories:
            fractions[category] = {}
            dm = rmssmfile.Get(category)
            for process in processes:
                if args.entries:
                    fractions[category][process] = dm.Get(process).GetEntries()
                else:
                    fractions[category][process] = dm.Get(process).Integral()
        output_data["_fraction"] = {}
        for process in processes:
            sm_signals = sum([output_data[a][process] for a in categories])
            mssm_signals = sum([fractions[a][process] for a in categories])
            output_data["_fraction"][process] = sm_signals / (sm_signals + mssm_signals)
    df = pd.DataFrame(output_data,
                      index=processes,
                      columns=categories + ["mass"] + ["_fraction"])
    return df, processes


def plot_contribution(data, channel, era, threshold, do_mask, entries):
    for signal_type in [
            "ggh_t", "ggh_b", "ggh_i", "ggH_t", "ggH_b", "ggH_i", "ggA_t",
            "ggA_b", "ggA_i", "bbH"
    ]:
        plt.figure(figsize=(9, 11))
        subset = data.filter(
            like=signal_type,
            axis=0).rename(index=lambda x: x.strip("{}_".format(signal_type)))
        subset.sort_values(by=['mass'], ascending=True, inplace=True)
        subset.drop(columns=["mass"], inplace=True)
        if do_mask:
            mask = subset.mask(subset < threshold,
                               True).mask(subset >= threshold, False)
            ax = sns.heatmap(subset,
                             cmap="OrRd",
                             annot=True,
                             linewidths=.5,
                             linecolor="grey",
                             annot_kws={'size': 12},
                             cbar_kws={"shrink": .8},
                             fmt=".3f",
                             mask=mask)
        else:
            ax = sns.heatmap(subset,
                             cmap="OrRd",
                             annot=True,
                             linewidths=.5,
                             linecolor="grey",
                             annot_kws={'size': 12},
                             cbar_kws={"shrink": .8},
                             fmt=".3f")
        xticks_labels = [
            catgegory.split("_")[1] for catgegory in subset.columns
        ]
        plt.xticks(np.arange(len(xticks_labels)) + .5, xticks_labels)
        if entries:
            plt.title("Signal Entries {}, {}, {}".format(era, channel,
                                                       signal_type),
                  loc='left',
                  fontsize=18)
        else:
            plt.title("Signal Integrals {}, {}, {}".format(era, channel,
                                                       signal_type),
                  loc='left',
                  fontsize=18)
        plt.ylabel("Mass Point")

        outputfile = "{}_{}_{}".format(era, channel, signal_type)
        print("Saving {}".format(outputfile))
        plt.tight_layout()
        plt.savefig("{}.pdf".format(outputfile))
        plt.savefig("{}.png".format(outputfile))
        plt.clf()


channel = args.channel
era = args.era
categories = sorted(
    ['{}_{}'.format(channel, cat) for cat in category_dict[channel]])

data, processes = get_integrals(channel, args.inputfile, args.mssmfile, categories, args.add_fraction, args.entries)
plot_contribution(data, channel, era, args.threshold, args.mask, args.entries)