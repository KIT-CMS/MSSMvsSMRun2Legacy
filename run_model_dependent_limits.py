import os 
import argparse 

# python run/model_dependent_limits.py --channel=fake --year=all
 
parser = argparse.ArgumentParser() 
parser.add_argument('--channel',help= 'Channel to run limits for', default='mt') 
parser.add_argument('--output','-o', help= 'Name of output directory', default='output') 
parser.add_argument('--year', help= 'Name of input year', default='2018') 
args = parser.parse_args() 
 
channel = args.channel 
output = args.output 
year = args.year 

analysis = 'mssm_vs_sm_classic'
 
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
   
 
variable="mt_tot_puppi" 

directory = "model_dependent_limits/%(output)s_%(analysis)s" % vars()


### Datacard creation ###
print "morph_parallel.py --output model_dependent_limits/%(output)s --analysis %(analysis)s --eras %(year_text)s --category_list input/%(cat_file)s --variable %(variable)s --parallel 5" % vars()
os.system("morph_parallel.py --output model_dependent_limits/%(output)s --analysis mssm_vs_sm_h125 --eras %(year_text)s --category_list input/%(cat_file)s --variable %(variable)s --parallel 5" % vars())
os.system("mkdir -p %(directory)s/combined/cmb/; rsync -av --progress %(directory)s/201?/htt_*/*  %(directory)s/combined/cmb/" % vars())

### Workspace creation ###
os.system("ulimit -s unlimited")
os.system("combineTool.py -M T2W -o ws_mh125.root  -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM --PO filePrefix=data/ --PO modelFile=13,Run2017,mh125_13.root --PO MSSM-NLO-Workspace=data/higgs_pt_v3_mssm_mode.root -i %(directory)s/combined/cmb/ --PO minTemplateMass=110.0 --PO maxTemplateMass=3200.0" % vars())

### Computing limits ###
os.system("mkdir -p %(directory)s/calculation_mh125_mssm_vs_sm_classic_h125")

os.system("ulimit -s unlimited")
os.system("combineTool.py -M AsymptoticGrid CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json -d %(directory)s/combined/cmb/ws_mh125.root --job-mode 'SGE' --prefix-file ic --sub-opts \"-q hep.q -l h_rt=3:0:0\" --task-name 'mssm_mh125_mssm_vs_sm_classic_h125_1' --redefineSignalPOI x --setParameters r=1 --freezeParameters r -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01" % vars())

