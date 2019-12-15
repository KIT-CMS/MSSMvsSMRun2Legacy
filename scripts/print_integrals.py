#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import ROOT as r
import sys
import os
import argparse
from prettytable import PrettyTable

parser = argparse.ArgumentParser(
    description="Compare Integrals of Processes between ML and Cutbased shapes"
)

parser.add_argument(
    "--inputfile",
    type=str,
    required=True,
    help="path to first shapes")
parser.add_argument(
    '--channel',
   type=str,
   required=True,
   help="Channel to analyze")

args = parser.parse_args()
r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

category_dict = {
    'cut': {
        'mt': [
            "wjets_control",

            "nobtag_boosted_0jet",
            "nobtag_0jet",

            "nobtag_1jet_lowpt_lowdeltar",
            "nobtag_1jet_lowpt_highdeltar",
            "nobtag_1jet_mediumpt_lowdeltar",
            "nobtag_1jet_mediumpt_highdeltar",
            "nobtag_1jet_highpt_lowdeltar",
            "nobtag_1jet_highpt_highdeltar",

            "nobtag_2jet_lowmjj_lowdeltar_lowjdeta",
            "nobtag_2jet_lowmjj_lowdeltar_highjdeta",
            "nobtag_2jet_lowmjj_highdeltar_lowjdeta",
            "nobtag_2jet_lowmjj_highdeltar_highjdeta",
            "nobtag_2jet_mediummjj_lowdeltar_lowjdeta",
            "nobtag_2jet_mediummjj_lowdeltar_highjdeta",
            "nobtag_2jet_mediummjj_highdeltar_lowjdeta",
            "nobtag_2jet_mediummjj_highdeltar_highjdeta",
            "nobtag_2jet_highmjj_lowdeltar",
            "nobtag_2jet_highmjj_highdeltar",

            "btag_tightmt",
            "btag_loosemt",
        ],
        'et': [
            "wjets_control",

            "nobtag_boosted_0jet",
            "nobtag_0jet",

            "nobtag_1jet_lowpt_lowdeltar",
            "nobtag_1jet_lowpt_highdeltar",
            "nobtag_1jet_mediumpt_lowdeltar",
            "nobtag_1jet_mediumpt_highdeltar",
            "nobtag_1jet_highpt_lowdeltar",
            "nobtag_1jet_highpt_highdeltar",

            "nobtag_2jet_lowmjj_lowdeltar_lowjdeta",
            "nobtag_2jet_lowmjj_lowdeltar_highjdeta",
            "nobtag_2jet_lowmjj_highdeltar_lowjdeta",
            "nobtag_2jet_lowmjj_highdeltar_highjdeta",
            "nobtag_2jet_mediummjj_lowdeltar_lowjdeta",
            "nobtag_2jet_mediummjj_lowdeltar_highjdeta",
            "nobtag_2jet_mediummjj_highdeltar_lowjdeta",
            "nobtag_2jet_mediummjj_highdeltar_highjdeta",
            "nobtag_2jet_highmjj_lowdeltar",
            "nobtag_2jet_highmjj_highdeltar",

            "btag_tightmt",
            "btag_loosemt",
        ],
        'tt': [
            "nobtag_boosted_0jet_lowdeltar",
            "nobtag_boosted_0jet_highdeltar",
            "nobtag_0jet_lowdeltar",
            "nobtag_0jet_highdeltar",

            "nobtag_1jet_lowpt_lowdeltar",
            "nobtag_1jet_lowpt_mediumdeltar",
            "nobtag_1jet_lowpt_highdeltar",
            "nobtag_1jet_mediumpt_lowdeltar",
            "nobtag_1jet_mediumpt_mediumdeltar",
            "nobtag_1jet_mediumpt_highdeltar",
            "nobtag_1jet_highpt_lowdeltar",
            "nobtag_1jet_highpt_mediumdeltar",
            "nobtag_1jet_highpt_highdeltar",

            "nobtag_2jet_lowmjj_lowdeltar_lowjdeta",
            "nobtag_2jet_lowmjj_lowdeltar_highjdeta",
            "nobtag_2jet_lowmjj_highdeltar_lowjdeta",
            "nobtag_2jet_lowmjj_highdeltar_highjdeta",
            "nobtag_2jet_mediummjj_lowdeltar_lowjdeta",
            "nobtag_2jet_mediummjj_lowdeltar_highjdeta",
            "nobtag_2jet_mediummjj_highdeltar_lowjdeta",
            "nobtag_2jet_mediummjj_highdeltar_highjdeta",
            "nobtag_2jet_highmjj_lowdeltar",
            "nobtag_2jet_highmjj_highdeltar",

            "btag",
        ],
        'em': [
            "ttbar_control",

            "nobtag_boosted_0jet",
            "nobtag_0jet",

            "nobtag_1jet_lowpt",
            "nobtag_1jet_mediumpt",
            "nobtag_1jet_highpt",

            "nobtag_2jet_lowmjj_lowjdeta",
            "nobtag_2jet_lowmjj_highjdeta",
            "nobtag_2jet_mediummjj_lowjdeta",
            "nobtag_2jet_mediummjj_highjdeta",
            "nobtag_2jet_highmjj",

            "btag_highdzeta",
            "btag_mediumdzeta",
            "btag_lowdzeta",
        ]
    }
}


def get_integrals(channel, f):
    categories = sorted([
        '{}_{}'.format(channel, cat) for cat in category_dict["cut"][channel]
    ])
    rfile = r.TFile.Open(f, "read")
    print("Looking into",f)
    processes = sorted([
        p.GetName() for p in rfile.Get(categories[0]).GetListOfKeys()
        if "Up" not in p.GetName() and "Down" not in p.GetName() and "data_obs" not in p.GetName() and "2H" not in p.GetName()
    ])
    for category in categories:
        print("\tCategory:",category)
        d = rfile.Get(category)
        for process in processes:
            integral = d.Get(process).Integral()
            print("\t\tProcess:",process,"Integral:",integral)
    rfile.Close()

channel = args.channel

get_integrals(channel, args.inputfile)
