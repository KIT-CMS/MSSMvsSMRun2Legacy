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
#    channel_perm = ['em','et','mt','tt','all']
    channel_perm = ['et','mt','tt']
  else:
    channel_perm = channel.split(",").append(channel)

  if year == 'all':
 #   year_perm = ['2016','2017','2018','all']
    year_perm = ['all']
  else:
    year_perm = year.split(',').append(year)
else:
  channel_perm = [channel]
  year_perm = [year]

  year_text = "2016,2017,2018"
  ### Datacard creation ###
  cat_file = 'mssm_classic_categories_2d_to_1d.txt'
  dc_creation_cmd = 'morph_parallel.py --output model_independent_limits/%(output)s_all_all --analysis "%(analysis)s" --eras %(year_text)s --category-list input/%(cat_file)s --variable "m_sv_VS_pt_tt_splitpT" --sm-gg-fractions data/higgs_pt_reweighting_fullRun2_v2.root --parallel 5 --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=1 " --sub-analysis "sm-like-light" --hSM-treatment "hSM-in-bg" --categorization="lowmass" --sm-like-hists="sm125" ' % vars()
  os.system(dc_creation_cmd)

  # take btag cats from m_sv binned one
  cat_file_btag = 'mssm_classic_categories_1d_btag.txt'
  dc_creation_cmd_2 = 'morph_parallel.py --output model_independent_limits/%(output)s_all_all --analysis "%(analysis)s" --eras %(year_text)s --category-list input/%(cat_file_btag)s --variable "m_sv_puppi" --sm-gg-fractions data/higgs_pt_reweighting_fullRun2_v2.root --parallel 5 --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=1 " --sub-analysis "sm-like-light" --hSM-treatment "hSM-in-bg" --categorization="lowmass" --sm-like-hists="sm125" ' % vars()
  os.system(dc_creation_cmd_2)


  cat_file_cr = 'mssm_classic_categories_cr.txt'
  dc_creation_cmd_3 = 'morph_parallel.py --output model_independent_limits/%(output)s_all_all --analysis "%(analysis)s" --eras %(year_text)s --category-list input/%(cat_file_cr)s --variable "mt_tot_puppi" --sm-gg-fractions data/higgs_pt_reweighting_fullRun2_v2.root --parallel 5 --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=1 " --sub-analysis "sm-like-light" --hSM-treatment "hSM-in-bg" --categorization="lowmass" --sm-like-hists="sm125" ' % vars()
  os.system(dc_creation_cmd_3)


for year in year_perm:
  for channel in channel_perm:
    ### Set up channel input ###
    
    ### Set up year input ###
    if year == "all":
      year_text = "2016,2017,2018"
      year_list = ["2016","2017","2018"]
    else:
      year_text = year
      year_list = year.split(",")

    #if year == channel: continue # skip full combination - temp, remove after!!
    #if not (year == 'all' or channel == 'all'): continue # don't run very fine breakdown by both channel and year for now 
    #if year == 'all': continue
    #if channel != 'mt': continue

    print year, channel
    
    directory = "model_independent_limits/%(output)s_%(channel)s_%(year)s_%(analysis)s" % vars()
    directory = "model_independent_limits/%(output)s_all_all_%(analysis)s" % vars()
    
    for yr in year_list:
      os.system("mkdir -p %(directory)s/%(yr)s/cmb/; rsync -av --progress %(directory)s/%(yr)s/htt_*/*  %(directory)s/%(yr)s/cmb/" % vars())
    os.system("mkdir -p %(directory)s/combined/cmb/; rsync -av --progress %(directory)s/201?/htt_*/*  %(directory)s/combined/cmb/" % vars())

    if channel !='all' and year !='all':
      os.system("mkdir -p %(directory)s/%(year)s/%(channel)s/; rsync -av --progress %(directory)s/%(year)s/htt_*/*  %(directory)s/%(year)s/%(channel)s/" % vars())
   
    if year == 'all' and channel != 'all':
      os.system("mkdir -p %(directory)s/combined/cmb/; rsync -av --progress %(directory)s/201?/htt_%(channel)s*/*  %(directory)s/combined/%(channel)s/; rsync -av --progress %(directory)s/201?/htt_em_2_*/*  %(directory)s/combined/%(channel)s/" % vars()) 
 
    ### Workspace creation ###
    
    os.system("ulimit -s unlimited")
   
    channel_str = channel
    year_str= year
    if channel == 'all': channel_str='cmb'
    if year == 'all': year_str = 'combined' 

    print year, channel, directory


    os.system('combineTool.py -M T2W -o "ws.root" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"'"'"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"'"'"' --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' --PO '"map=^.*/qqX$:r_qqX[0]"' --PO '"'"'"map=^.*/ggX_(i|t|b).?$:r_ggX[0,0,200]"'"'"'  -i %(directory)s/%(year_str)s/%(channel_str)s -m 95 --parallel 8' % vars())


    os.system('combineTool.py -m "60,65,70,75,80,85,90,95,100,105,110,120,125,130,140,160,180,200,250" -M Significance --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".ggH.v2" --task-name ggH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s_pvalue --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0"' % vars())

    os.system('combineTool.py -m "60,65,70,75,80,85,90,95,100,105,110,120,125,130,140,160,180,200,250" -M Significance --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_bbH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".bbH.v2" --task-name bbH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s_pvalue --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0"' % vars())

    os.system('combineTool.py -m "60,80,95,100,120,125,130,140,160,180,200,250" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".ggH.v2" --task-name ggH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" ' % vars())

    os.system('combineTool.py -m "60,80,95,100,120,125,130,140,160,180,200,250" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_bbH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".bbH.v2" --task-name bbH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" ' % vars())

#    # azimov fits:
#    
#    os.system('combineTool.py -m "100" -M GenerateOnly --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=5.8,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".ggH.Toys100"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01  -t -1 --saveToys' % vars())
#
#    os.system('combineTool.py -m "60,65,70,75,80,85,90,95,100,105,110,120,125,130,140,160,180,200,250" -M Significance --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".ggH.azimov" --task-name ggH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s_pvalue --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t -1 --toysFile /vols/cms/dw515/MSSMLowMass/CMSSW_10_2_25/src/CombineHarvester/MSSMvsSMRun2Legacy/%(directory)s/%(year_str)s/%(channel_str)s/higgsCombine.ggH.Toys100.GenerateOnly.mH100.123456.root ' % vars())
#
#    os.system('combineTool.py -m "60,65,70,75,80,85,90,95,100,105,110,120,125,130,140,160,180,200,250" -M Significance --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_bbH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".bbH.azimov" --task-name bbH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s_pvalue --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t -1 --toysFile /vols/cms/dw515/MSSMLowMass/CMSSW_10_2_25/src/CombineHarvester/MSSMvsSMRun2Legacy/%(directory)s/%(year_str)s/%(channel_str)s/higgsCombine.ggH.Toys100.GenerateOnly.mH100.123456.root ' % vars())
#
#    os.system('combineTool.py -m "60,80,95,100,120,125,130,140,160,180,200,250" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".ggH.azimov" --task-name ggH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t -1 --toysFile /vols/cms/dw515/MSSMLowMass/CMSSW_10_2_25/src/CombineHarvester/MSSMvsSMRun2Legacy/%(directory)s/%(year_str)s/%(channel_str)s/higgsCombine.ggH.Toys100.GenerateOnly.mH100.123456.root ' % vars())
#
#    os.system('combineTool.py -m "60,80,95,100,120,125,130,140,160,180,200,250" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_bbH -d %(directory)s/%(year_str)s/%(channel_str)s/ws.root --there -n ".bbH.azimov" --task-name bbH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t -1 --toysFile /vols/cms/dw515/MSSMLowMass/CMSSW_10_2_25/src/CombineHarvester/MSSMvsSMRun2Legacy/%(directory)s/%(year_str)s/%(channel_str)s/higgsCombine.ggH.Toys100.GenerateOnly.mH100.123456.root ' % vars())
