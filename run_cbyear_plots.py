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
    

for year in year_perm:
  for channel in channel_perm:
    ### Set up channel input ###
    if channel == "all":
      cat_file = 'mssm_classic_categories_nobtag.txt'
    elif channel == "fake":
      cat_file = 'mssm_fake_categories.txt'
    elif channel in ["et","tt","mt","em"]:
      cat_file = 'mssm_%(channel)s_categories_nobtag.txt' % vars()
    
    ### Set up year input ###
    if year == "all":
      year_text = "2016,2017,2018"
      year_list = ["2016","2017","2018"]
    else:
      year_text = year
      year_list = year.split(",")
    
    ### Datacard creation ###
    dc_creation_cmd = 'morph_parallel.py --output model_independent_limits/%(output)s_%(channel)s_%(year)s --analysis "%(analysis)s" --eras %(year_text)s --category-list input/%(cat_file)s --variable "mt_tot_puppi" --sm-gg-fractions data/higgs_pt_reweighting_fullRun2_v2.root --parallel 5 --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=1 --cbyear_plot=true " --sub-analysis "sm-like-light" --hSM-treatment "hSM-in-bg" --categorization="classic" --sm-like-hists="sm125" ' % vars()
    os.system(dc_creation_cmd)
    # take btag cats from m_sv binned one
    cat_file_btag = 'mssm_classic_categories_btag.txt'
    dc_creation_cmd_2 = 'morph_parallel.py --output model_independent_limits/%(output)s_%(channel)s_%(year)s --analysis "%(analysis)s" --eras %(year_text)s --category-list input/%(cat_file_btag)s --variable "mt_tot_puppi" --sm-gg-fractions data/higgs_pt_reweighting_fullRun2_v2.root --parallel 5 --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=1 --cbyear_plot=true " --sub-analysis "sm-like-light" --hSM-treatment "hSM-in-bg" --categorization="classic" --sm-like-hists="sm125" ' % vars()
    os.system(dc_creation_cmd_2)

    cat_file = 'mssm_classic_categories_2d_to_1d.txt'
    dc_creation_cmd = 'morph_parallel.py --output model_independent_limits/%(output)s_all_all --analysis "%(analysis)s" --eras %(year_text)s --category-list input/%(cat_file)s --variable "m_sv_VS_pt_tt_splitpT" --sm-gg-fractions data/higgs_pt_reweighting_fullRun2_v2.root --parallel 5 --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=1 --cbyear_plot=true " --sub-analysis "sm-like-light" --hSM-treatment "hSM-in-bg" --categorization="lowmass" --sm-like-hists="sm125" ' % vars()
    os.system(dc_creation_cmd)
    
    directory = "model_independent_limits/%(output)s_%(channel)s_%(year)s_%(analysis)s" % vars()

    for c in ['em','lt','tt']:
      for b in [32,35,432,332]:
        if b in [32,35]: 
          mass = 1200
          fit_directory="model_independent_limits/Jan12_mt_tot_all_all_bsm-model-indep/"
          freeze='--freeze r_bbH=0.00306'
        else: 
          mass = 100
          fit_directory="model_independent_limits/Feb10_all_all_bsm-model-indep/"
          freeze='--freeze r_bbH=5.758'

        print 'Doing bin: htt_%(c)s_%(b)s' % vars()
        os.system('combineTool.py -M T2W -o "ws.root" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"'"'"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"'"'"' --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' --PO '"map=^.*/qqX$:r_qqX[0]"' --PO '"'"'"map=^.*/ggX_(i|t|b).?$:r_ggX[0,0,200]"'"'"'  -i %(directory)s/htt_%(c)s_%(b)s -m %(mass)s --parallel 8' % vars())
        os.system('PostFitShapesFromWorkspace -w %(directory)s/htt_%(c)s_%(b)s/ws.root -d %(directory)s/htt_%(c)s_%(b)s/combined.txt.cmb --fitresult %(fit_directory)s/combined/cmb/multidimfitggH.m%(mass)s.bestfit.robustHesse.root:fit_mdf -o shapes_cbyears_%(c)s_%(b)s.root --skip-prefit=true  --mass %(mass)s --total-shapes=true --postfit %(freeze)s' % vars())

