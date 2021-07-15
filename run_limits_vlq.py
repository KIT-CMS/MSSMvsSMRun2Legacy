import os
import argparse

# python run_model_independent_limits.py --channel=fake --year=2016 

parser = argparse.ArgumentParser()
parser.add_argument('--channel',help= 'Channel to run limits for', default='mt')
parser.add_argument('--output','-o', help= 'Name of output directory', default='output')
parser.add_argument('--year', help= 'Name of input year', default='2018')
parser.add_argument('--combine', help= 'Name of input year', default=False)
parser.add_argument('--all_perm', help= 'Run all permutations of year and channel inputs', default=False)
args = parser.parse_args()

channel = args.channel
output = args.output
year = args.year
all_perm = args.all_perm
combine = args.combine

def CreateBatchJob(name,cmd_list):
  if os.path.exists(name): os.system('rm %(name)s' % vars())
  os.system('echo "#!/bin/bash" >> %(name)s' % vars())
  os.system('echo "ulimit -s unlimited" >> %(name)s' % vars())
  os.system('echo "cd /vols/cms/gu18/CH/CMSSW_10_2_25/src" >> %(name)s' % vars())
  os.system('echo "export SCRAM_ARCH=slc7_amd64_gcc700" >> %(name)s' % vars())
  os.system('echo "source /vols/grid/cms/setup.sh" >> %(name)s' % vars())
  os.system('echo "eval \'scramv1 runtime -sh\'" >> %(name)s' % vars())
  os.system('echo "cd /vols/cms/gu18/CH/CMSSW_10_2_25/src/CombineHarvester/MSSMvsSMRun2Legacy" >> %(name)s' % vars())
  for cmd in cmd_list:
    os.system('echo "%(cmd)s" >> %(name)s' % vars())
  os.system('chmod +x %(name)s' % vars())
  print "Created job:",name

def SubmitBatchJob(name,time=180,memory=24,cores=1):
  error_log = name.replace('.sh','_error.log')
  output_log = name.replace('.sh','_output.log')
  if os.path.exists(error_log): os.system('rm %(error_log)s' % vars())
  if os.path.exists(output_log): os.system('rm %(output_log)s' % vars())
  if cores>1: os.system('qsub -e %(error_log)s -o %(output_log)s -V -q hep.q -pe hep.pe %(cores)s -l h_rt=0:%(time)s:0 -l h_vmem=%(memory)sG -cwd %(name)s' % vars())
  else: os.system('qsub -e %(error_log)s -o %(output_log)s -V -q hep.q -l h_rt=0:%(time)s:0 -l h_vmem=%(memory)sG -cwd %(name)s' % vars())

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
      
    
    cmds = []

    ### Datacard creation ###
    dc_creation_cmd = 'morph_parallel.py --output model_independent_limits/%(output)s_%(channel)s_%(year)s --analysis "%(analysis)s" --eras %(year_text)s --category_list input/%(cat_file)s --variable "mt_tot_puppi" --parallel 5 --additional_arguments="--auto_rebin=1"' % vars()
    cmds.append(dc_creation_cmd)

    
    directory = "model_independent_limits/%(output)s_%(channel)s_%(year)s_%(analysis)s" % vars()
    
    for yr in year_list:
      cmds.append("mkdir -p %(directory)s/%(yr)s/cmb/; rsync -av --progress %(directory)s/%(yr)s/htt_*/*  %(directory)s/%(yr)s/cmb/" % vars()) 
    cmds.append("mkdir -p %(directory)s/combined/cmb/; rsync -av --progress %(directory)s/201?/htt_*/*  %(directory)s/combined/cmb/" % vars())
    
    ### Workspace creation ###
    
    cmds.append("ulimit -s unlimited")
    cmds.append('combineTool.py -M T2W -o \\"ws.root\\" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO \'\\"map=^.*/VLQ_betaRd33_0_M-.?$:r_vlq[0,0,250]\\"\' -i %(directory)s/{%(year_text)s,combined}/cmb -m 2000 --parallel 4' % vars())

    CreateBatchJob("model_indep.sh",cmds)
    SubmitBatchJob("model_indep.sh")

    ### Run model-independent limits ###

    if combine:

      if year == "all" and channel == "all": time = '6'
      else: time ='3'
 
      os.system('combineTool.py -m "2000,3000,4000" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist input/vlq_boundaries.json -d %(directory)s/combined/cmb/ws.root --there -n ".VLQ_betaRd33_0_M-" --task-name ggH_full_combined_%(analysis)s_%(channel)s_%(year)s_%(output)s --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1 --job-mode \'SGE\' --prefix-file ic --sub-opts "-q hep.q -l h_rt=%(time)s:0:0"' % vars())

    if not combine:
      CreateBatchJob("model_indep_%(channel)s_%(year)s.sh" % vars(),cmds)
      SubmitBatchJob("model_indep_%(channel)s_%(year)s.sh" % vars())
    

