#! /usr/bin/env python3

import ROOT as r
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Script to determine outliers in Graphs of MSSM benchmark CLs results by checking a linear interpolation of neighbour points.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--input", required=True, help="input ROOT file including the CLs results, usually asymptotic_grid.root")
parser.add_argument("--output", required=True, help="output ROOT file including the CLs results corrected for outliers")
parser.add_argument("--outlier-threshold", type=float, default=3.0, help="Threshold factor how much the outlier is allowed to deviate from linearly extrapolated value.")
parser.add_argument("--continuity-threshold", type=float, default=1.5, help="Factor, how much neighbour points are allow to deviate from their mean.")
parser.add_argument("--consideration-threshold", type=float, default=0.001, help="CLs value from which on one should consider numerical values to avoid numerical fluctuations.")

args = parser.parse_args()

f = r.TFile.Open(args.input, "read")

levels = {
    "exp-2" : None,
    "exp-1" : None,
    "exp0" : None,
    "exp+1" : None,
    "exp+2" : None,
    "obs" : None
}

for l in levels:
    graph = f.Get(l)
    npoints = graph.GetN()

    x = graph.GetX()
    y = graph.GetY()
    z = graph.GetZ()
    d = { "mA" : [x[i] for i in range(npoints)], "tanb": [y[i] for i in range(npoints)], "CLs" : [z[i] for i in range(npoints)]}
    levels[l] = pd.DataFrame(data=d)
    print(f"Level {l} with {len(levels[l].index)} points")

f.Close()

problematic_points = set()

for level_name, graph in levels.items():
    print(f"Showing {level_name}, scanning mA:")
    tanb_vals = sorted(set(graph["tanb"]))
    mA_vals = sorted(set(graph["mA"]))
    for mA in mA_vals:
        graph_mA = graph[graph["mA"] == mA]
        graph_tanb_vals = sorted(graph_mA["tanb"].values)
        n_tanb_triplets = len(graph_tanb_vals) - 2
        for i in range(n_tanb_triplets):
            tanb_left = graph_tanb_vals[i]
            tanb_center = graph_tanb_vals[i+1]
            tanb_right = graph_tanb_vals[i+2]
            CLs_left = graph_mA[graph_mA["tanb"] == tanb_left]["CLs"].values[0]
            CLs_center = graph_mA[graph_mA["tanb"] == tanb_center]["CLs"].values[0]
            CLs_right = graph_mA[graph_mA["tanb"] == tanb_right]["CLs"].values[0]
            CLs_linear_extrapolated = (CLs_right - CLs_left)/(tanb_right - tanb_left) * (tanb_center - tanb_left) + CLs_left
            bigger_outlier_check = abs(CLs_center / CLs_linear_extrapolated) > args.outlier_threshold and any(np.array([CLs_left, CLs_center, CLs_right]) >= args.consideration_threshold)
            smaller_outlier_check = abs(CLs_center / CLs_linear_extrapolated) < 1/args.outlier_threshold and any(np.array([CLs_left, CLs_center, CLs_right]) >= args.consideration_threshold)
            continuity_check = abs(2 * CLs_left / (CLs_left + CLs_right)) < args.continuity_threshold and abs(2 * CLs_left / (CLs_left + CLs_right)) > 1/args.continuity_threshold
            if bigger_outlier_check and continuity_check:
                print(f"\tOutlier: mA = {mA}, tanb = {tanb_center}, CLs = {CLs_center}; neighbour values: CLs left = {CLs_left}, CLs right = {CLs_right}")
                problematic_points.add((mA, tanb_center))
            elif smaller_outlier_check and continuity_check:
                print(f"\tOutlier: mA = {mA}, tanb = {tanb_center}, CLs = {CLs_center}; neighbour values: CLs left = {CLs_left}, CLs right = {CLs_right}")
                problematic_points.add((mA, tanb_center))
    print(f"Showing {level_name}, scanning tanb:")
    for tanb in tanb_vals:
        graph_tanb = graph[graph["tanb"] == tanb]
        graph_mA_vals = sorted(graph_tanb["mA"].values)
        n_mA_triplets = len(graph_mA_vals) - 2
        for i in range(n_mA_triplets):
            mA_left = graph_mA_vals[i]
            mA_center = graph_mA_vals[i+1]
            mA_right = graph_mA_vals[i+2]
            CLs_left = graph_tanb[graph_tanb["mA"] == mA_left]["CLs"].values[0]
            CLs_center = graph_tanb[graph_tanb["mA"] == mA_center]["CLs"].values[0]
            CLs_right = graph_tanb[graph_tanb["mA"] == mA_right]["CLs"].values[0]
            CLs_linear_extrapolated = (CLs_right - CLs_left)/(tanb_right - tanb_left) * (tanb_center - tanb_left) + CLs_left
            bigger_outlier_check = abs(CLs_center / CLs_linear_extrapolated) > args.outlier_threshold and any(np.array([CLs_left, CLs_center, CLs_right]) >= args.consideration_threshold)
            smaller_outlier_check = abs(CLs_center / CLs_linear_extrapolated) < 1/args.outlier_threshold and any(np.array([CLs_left, CLs_center, CLs_right]) >= args.consideration_threshold)
            continuity_check = abs(2 * CLs_left / (CLs_left + CLs_right)) < args.continuity_threshold and abs(2 * CLs_left / (CLs_left + CLs_right)) > 1/args.continuity_threshold
            if bigger_outlier_check and continuity_check:
                print(f"\tOutlier: mA = {mA_center}, tanb = {tanb}, CLs = {CLs_center}; neighbour values: CLs left = {CLs_left}, CLs right = {CLs_right}")
                problematic_points.add((mA_center, tanb))
            elif smaller_outlier_check and continuity_check:
                print(f"\tOutlier: mA = {mA_center}, tanb = {tanb}, CLs = {CLs_center}; neighbour values: CLs left = {CLs_left}, CLs right = {CLs_right}")
                problematic_points.add((mA_center, tanb))

problematic_points = sorted(problematic_points)
print(f"Problematic points: {problematic_points}")
print(f"In total: {len(problematic_points)}")

for level_name in levels:
    for p in problematic_points:
        point_filter = levels[level_name]["mA"] == p[0]
        point_filter *= levels[level_name]["tanb"] == p[1]
        point_filter = np.logical_not(point_filter)
        levels[level_name] = levels[level_name][point_filter]

out = r.TFile.Open(args.output, "recreate")
out.cd()

for level_name, graph in levels.items():
    n_points = len(graph.index)
    print(f"Level {level_name} reduced to {n_points} points")
    new_graph = r.TGraph2D(level_name, level_name, n_points, graph["mA"].values, graph["tanb"].values, graph["CLs"].values)
    new_graph.Write()
