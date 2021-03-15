import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--channel',help= 'Channel to run limits for', default='mt')
parser.add_argument('--output','-o', help= 'Name of output directory', default='output')
parser.add_argument('--year', help= 'Name of input year', default='2018')
parser.add_argument('--all_perm', help= 'Run all permutations of year and channel inputs', default=False)
args = parser.parse_args()

channel = args.channel
output = args.output
year = args.year
all_perm = args.all_perm

analysis = "mssm_classic"

if all_perm:
  if channel == 'fake':
    channel_perm = ['et','mt','tt','fake']
  elif channel == 'all':
    channel_perm = ['em','et','mt','tt','all']
  else:
    channel_perm = channel.split(",").append(channel)

  if year == 'all':
    year_perm = ['2016','2017','2018','all']
  else:
    year_perm = year.split(',').append(year)
else:
  channel_perm = [channel]
  year_perm = [year]

for year in year_perm:
  for channel in channel_perm:

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
    
    
    directory = "model_independent_limits/%(output)s_%(channel)s_%(year)s_%(analysis)s" % vars()
    
    
    ### Collecting limits ###
    for p in ["gg","bb"]:
      print "combineTool.py -M CollectLimits %(directory)s/combined/cmb/higgsCombine.%(p)sH*.root --use-dirs -o %(directory)s/combined/cmb/mssm_%(p)sH_combined.json" % vars()
      os.system("combineTool.py -M CollectLimits %(directory)s/combined/cmb/higgsCombine.%(p)sH*.root --use-dirs -o %(directory)s/combined/cmb/mssm_%(p)sH_combined.json" % vars())
    
    ### Plotting limits ###
    for p in ["gg","bb"]:
      os.system('plotMSSMLimits.py --cms-sub "Preliminary" --title-right "137 fb^{-1} (13 TeV)" --process "%(p)s#phi" --y-axis-min 0.0001 --y-axis-max 1000.0 --show obs,exp %(directory)s/combined/cmb/mssm_%(p)sH_combined_cmb.json  --output %(directory)s/model-independent_combined_%(p)sH_cmb --logx --logy' % vars())
    

