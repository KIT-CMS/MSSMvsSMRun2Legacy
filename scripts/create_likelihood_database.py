#!/usr/bin/env python

import argparse
import itertools
import operator
import logging
from array import array

import pandas as pd

import ROOT
ROOT.gROOT.SetBatch()
ROOT.PyConfig.IgnoreCommandLineOptions = True

import CombineHarvester.CombineTools.plotting as plot


logger = logging.getLogger("")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("files",
                        nargs="+",
                        help="Input files")
    parser.add_argument("--output", "-o",
                        default="limit",
                        help="Name of the output plot without file extension")
    return parser.parse_args()


def setup_logging(level=logging.INFO):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def rezero_tgraph2d(graph, perform=True):
    # Mostly copied from plot.ReZeroTGraph
    fit_x = 0.
    fit_y = 0.
    fit_z = 0.
    # Find minimum of fit as first value with
    for i in xrange(graph.GetN()):
        if graph.GetZ()[i] == 0.:
            fit_x = graph.GetX()[i]
            fit_y = graph.GetY()[i]
            fit_z = graph.GetZ()[i]
            break
    min_x = 0.
    min_y = 0.
    min_z = 0.
    for i in xrange(graph.GetN()):
        if graph.GetZ()[i] < min_z:
            min_z = graph.GetZ()[i]
            min_y = graph.GetY()[i]
            min_x = graph.GetX()[i]
    if min_z < fit_z:
        logging.info('[ReZeroTGraph] Fit minimum was (%f, %f, %f)' % (fit_x, fit_y, fit_z))
        logging.info('[ReZeroTGraph] Better minimum was (%f, %f, %f)' % (min_x, min_y, min_z))
        if perform:
            for i in xrange(graph.GetN()):
                before = graph.GetZ()[i]
                graph.GetZ()[i] -= min_z
                after = graph.GetZ()[i]
                # print 'Point %i, before=%f, after=%f' % (i, before, after)
    return min_z


def convert_graph_to_dataframe(graph):
    # Get all x and y values from graph
    logger.info("Starting conversion from TGraph to dataframe...")
    x_vals = set()
    y_vals = set()
    points = []
    gx = graph.GetX()
    gy = graph.GetY()
    gz = graph.GetZ()
    logger.debug("Retrieving information on scan points from TGraph...")
    for i in xrange(graph.GetN()):
        x_vals.add(gx[i])
        y_vals.add(gy[i])
        points.append((gx[i], gy[i], gz[i]))

    x_vals = sorted(x_vals)
    y_vals = sorted(y_vals)
    points = sorted(points, key=operator.itemgetter(0,1))

    # For values not contained in graph, interpolate the values
    # using the interpolation implemented in TGraph
    logger.debug("Searching for entries missing in the created TGraph...")
    getxy = operator.itemgetter(0,1)
    existing_points = map(getxy, points)
    missing = set(itertools.product(x_vals, y_vals)) - set(existing_points)
    logger.info("Found {} missing entries in scan".format(len(missing)))
    if len(missing) > 0:
        logger.info("Will set their values to the interpolated ones...")
    miss_entries = map(lambda x: (x[0], x[1], graph.Interpolate(*x)), missing)

    # Build dataframe from the uniqe x and y values
    df = pd.DataFrame(data=points,
                      columns=["r_ggH", "r_bbH", "deltaNLL"])
    logger.debug("Dataframe successfully created...")
    # Fill dataframe with values contained in graph
    df = df.append(pd.DataFrame(data=miss_entries,
                                columns=["r_ggH", "r_bbH", "deltaNLL"]),
                                ignore_index=True,
                                sort=True)
    logger.debug("Added missing values from scan to dataframe...")

    logger.debug("Sorting dataframe to restore correct order of scan...")
    df.sort_values(by=["r_ggH", "r_bbH"], inplace=True)
    return df


def convert_dataframe_to_tree(df, best_fit, offset=0.):
    # Create tree structure and branch pointers
    tree = ROOT.TTree("limit", "limit")
    ggh = array("f", [0.])
    tree.Branch("r_ggH", ggh, "r_ggH/F")
    bbh = array("f", [0.])
    tree.Branch("r_bbH", bbh, "r_bbH/F")
    deltaNLL = array("f", [0.])
    tree.Branch("deltaNLL", deltaNLL, "deltaNLL/F")
    quantileExp = array("f", [0.])
    tree.Branch("quantileExpected", quantileExp, "quantileExpected/F")
    # First get best fit point and write it to file
    x, y = ROOT.Double(), ROOT.Double()
    best_fit.GetPoint(0, x, y)
    ggh[0] = x
    bbh[0] = y
    deltaNLL[0] = 0. - offset
    quantileExp[0] = -1
    tree.Fill()
    # Loop over dataframe entries and fill the tree
    for index, row in df.iterrows():
        ggh[0] = row["r_ggH"]
        bbh[0] = row["r_bbH"]
        deltaNLL[0] = row["deltaNLL"]
        quantileExp[0] = ROOT.Math.chisquared_cdf_c(2*row["deltaNLL"], 2)
        tree.Fill()
    return tree


def main(args):
    # Get tree with scan values from all input files
    limit = plot.MakeTChain(args.files, 'limit')

    # Write Tree entries in TGraph to get structure for likelihood database
    graph = plot.TGraph2DFromTree(limit,
                                  "r_ggH", "r_bbH",
                                  'deltaNLL',
                                  'quantileExpected > -0.5 && deltaNLL < 1000')
                                  # 'quantileExpected > -0.5 && deltaNLL > 0 && deltaNLL < 1000')

    # rezero the TGraph to have sensible values in the output file
    min_delta_nll = rezero_tgraph2d(graph)
    df = convert_graph_to_dataframe(graph)
    df.to_csv(args.output, sep=" ",
              float_format="%.6f",
              header=False, index=False,
              na_rep="NaN",
              columns=["r_ggH", "r_bbH", "deltaNLL"])
    best = plot.TGraphFromTree(
        limit, "r_ggH", "r_bbH", 'deltaNLL == 0')
    plot.RemoveGraphXDuplicates(best)
    # Write back a root TTree with the modified contents.
    outfile = ROOT.TFile(args.output.replace(".txt", ".root"), "recreate")
    tree = convert_dataframe_to_tree(df, best, offset=min_delta_nll)
    # Add the best fit values from the limit trees
    outfile.Write()
    return


if __name__ == "__main__":
    args = parse_args()
    setup_logging(level=logging.INFO)
    main(args)
