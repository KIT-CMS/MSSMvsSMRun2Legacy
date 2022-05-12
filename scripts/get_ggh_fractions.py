#! /usr/bin/env python
import ROOT as r
import sys

workspace = sys.argv[1]
masses = [
  60., 80., 95., 100., 120., 125., 130., 140., 160., 180., 200.,
  250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000.,
  1200., 1400., 1600., 1800., 2000., 2300., 2600., 2900., 3200., 3500.
]

f = r.TFile.Open(workspace, "read")
ws = f.Get("w")

mass = ws.var("MH")
ggh_i_frac = ws.function("ggh_i_frac")
ggh_t_frac = ws.function("ggh_t_frac")
ggh_b_frac = ws.function("ggh_b_frac")

content = ["MH,ggh_t_frac,ggh_b_frac,ggh_i_frac\n"]

for m  in masses:
  mass.setVal(m)
  content.append(",".join([str(m),str(ggh_t_frac.evaluate()),str(ggh_b_frac.evaluate()),str(ggh_i_frac.evaluate())])+"\n")

with open("sm_ggh_fractions.csv", "w") as out:
  out.writelines(content)
