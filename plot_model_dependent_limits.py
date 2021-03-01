import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--channel',help= 'Channel to run limits for', default='mt')
parser.add_argument('--output','-o', help= 'Name of output directory', default='output')
parser.add_argument('--year', help= 'Name of input year', default='2018')
args = parser.parse_args()

channel = args.channel
output = args.output
year = args.year

### Set up channel input ### 
if channel == "all":
  cat_file = 'mssm_classic_categories.txt'
elif channel == "fake":
  cat_file = 'mssm_fake_categories.txt'
elif channel in ["et","tt","mt","em"]:
  cat_file = 'mssm_%(channel)s_classic_categories.txt' % vars()

### Set up year input ### 
if year == "all":
  year_text = "2016,2017,2018"
  year_list = ["2016","2017","2018"]
else:
  year_text = year
  year_list = year.split(",")


variable="mttot"

directory = "model_dependent_limits/%(output)s_mssm_vs_sm_classic_h125" % vars()

os.system("ulimit -s unlimited")

os.system("plotLimitGrid.py asymptotic_grid.root --scenario-label=\"M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)\" --output mssm_mh125_mssm_vs_sm_classic_h125  --title-right=\"137 fb^{-1} (13 TeV)\" --cms-sub=\"Preliminary\" --contours=\"exp-2,exp-1,exp0,exp+1,exp+2,obs\" --model_file=data/mh125_13.root --y-range 2.0,60.0 --x-title \"m_{A} [GeV]\"")
