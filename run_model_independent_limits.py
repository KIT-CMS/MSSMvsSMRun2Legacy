import os
import argparse

# python run_model_independent_limits.py --channel=fake --year=2016 

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

analysis = 'mssm_classic'

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
      cat_file = 'mssm_%(channel)s_categories.txt' % vars()
    
    ### Set up year input ###
    if year == "all":
      year_text = "2016,2017,2018"
      year_list = ["2016","2017","2018"]
    else:
      year_text = year
      year_list = year.split(",")
      
    
    ### Datacard creation ###
    dc_creation_cmd = 'morph_parallel.py --output model_independent_limits/%(output)s_%(channel)s_%(year)s --analysis "%(analysis)s" --eras %(year_text)s --category_list input/%(cat_file)s --variable "mt_tot_puppi" --sm_gg_fractions data/higgs_pt_reweighting_fullRun2.root --parallel 5 --additional_arguments="--auto_rebin=1"' % vars()
    print dc_creation_cmd
    os.system(dc_creation_cmd)
    
    directory = "model_independent_limits/%(output)s_%(channel)s_%(year)s_%(analysis)s" % vars()
    
    for yr in year_list:
      os.system("mkdir -p %(directory)s/%(yr)s/cmb/; rsync -av --progress %(directory)s/%(yr)s/htt_*/*  %(directory)s/%(yr)s/cmb/" % vars())
    os.system("mkdir -p %(directory)s/combined/cmb/; rsync -av --progress %(directory)s/201?/htt_*/*  %(directory)s/combined/cmb/" % vars())
    
    ### Workspace creation ###
    
    os.system("ulimit -s unlimited")
    os.system('combineTool.py -M T2W -o "ws.root" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"'"'"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"'"'"' --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' -i %(directory)s/{%(year_text)s,combined}/cmb -m 110 --parallel 4' % vars())
    
    ### Run model-independent limits ###

    os.system('combineTool.py -m "60,80,100,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist input/mssm_boundaries.json --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_ggH -d %(directory)s/combined/cmb/ws.root --there -n ".ggH" --task-name ggH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0"' % vars())

    os.system('combineTool.py -m "60,80,100,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist input/mssm_boundaries.json --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_bbH -d %(directory)s/combined/cmb/ws.root --there -n ".bbH" --task-name bbH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0"' % vars())
    

