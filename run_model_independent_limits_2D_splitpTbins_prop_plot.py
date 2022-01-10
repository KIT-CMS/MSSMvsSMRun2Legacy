import os
import argparse

# python run_model_independent_limits_unblinding.py -o unblinding 

parser = argparse.ArgumentParser()
parser.add_argument('--channel',help= 'Channel to run limits for', default='all')
parser.add_argument('--output','-o', help= 'Name of output directory', default='output')
parser.add_argument('--year', help= 'Name of input year', default='all')
parser.add_argument('--all_perm', help= 'Run all permutations of year and channel inputs', default=False)
args = parser.parse_args()

channel = args.channel
output = args.output
year = args.year
all_perm = args.all_perm

analysis = 'bsm-model-indep'

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

  year_text = "2016,2017,2018"
  ### Datacard creation ###
  cat_file = 'mssm_classic_categories_2d_to_1d.txt'
  dc_creation_cmd = 'morph_parallel.py --output model_independent_limits/%(output)s_all_all --analysis "%(analysis)s" --eras %(year_text)s --category-list input/%(cat_file)s --variable "m_sv_VS_pt_tt_splitpT" --sm-gg-fractions data/higgs_pt_reweighting_fullRun2_v2.root --parallel 5 --additional-arguments="--prop_plot=true --auto_rebin=1 --manual_rebin=1 --real_data=1 " --sub-analysis "sm-like-light" --hSM-treatment "hSM-in-bg" --categorization="lowmass" --sm-like-hists="sm125" ' % vars()


  os.system(dc_creation_cmd)
    

year='all'
channel='all'

for x in ["plot_cmb","plot_50to100","plot_100to200","plot_0to50","plot_GT200"]:

    ### Set up channel input ###
    
    directory = "model_independent_limits/%(output)s_all_all_%(analysis)s" % vars()
    
    os.system("mkdir -p %(directory)s/combined/%(x)s/; rsync -av --progress %(directory)s/201?/%(x)s/*  %(directory)s/combined/%(x)s/" % vars())
 
    ### Workspace creation ###
    
    os.system("ulimit -s unlimited")
   
    channel_str='%(x)s' % vars()
    year_str = 'combined' 

    os.system('combineTool.py -M T2W -o "ws.root" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"'"'"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"'"'"' --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' -i %(directory)s/%(year_str)s/%(channel_str)s -m 100 --parallel 8' % vars())

    # s+b fit
    os.system('PostFitShapesFromWorkspace -w %(directory)s/%(year_str)s/%(channel_str)s/ws.root -d %(directory)s/%(year_str)s/%(channel_str)s/combined.txt.cmb --fitresult %(directory)s/combined/cmb/multidimfitggH.m100.bestfit.robustHesse.root:fit_mdf -o shapes_prop_plot_postfit_${x}_dec20.root --skip-prefit=true   --mass 100 --total-shapes=true --postfit')

    # b-only fit
    os.system('PostFitShapesFromWorkspace -w %(directory)s/%(year_str)s/%(channel_str)s/ws.root -d %(directory)s/%(year_str)s/%(channel_str)s/combined.txt.cmb --fitresult %(directory)s/combined/cmb/multidimfitggH.bkgOnly.bestfit.robustHesse.root:fit_mdf -o shapes_prop_plot_postfit_bkgonly_${x}_dec20.root --skip-prefit=true   --mass 100 --total-shapes=true --postfit')


