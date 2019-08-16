#!/bin/bash

### Setup of CMSSW release
NUM_CORES=10
CMSSW=CMSSW_10_2_16_UL

export SCRAM_ARCH=slc7_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project $CMSSW; pushd $CMSSW/src
eval `scramv1 runtime -sh`

# combine on 102X slc7
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.0.1
cd -

# CombineHarvester (current master for 102X)
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester

# MSSM vs SM analysis specific code
git clone https://github.com/KIT-CMS/MSSMvsSMRun2Legacy.git CombineHarvester/MSSMvsSMRun2Legacy

# Install LHCHXSWGMSSMNeutral packages
wget -P CombineHarvester/MSSMvsSMRun2Legacy/python -O mssm_xs_tools.py https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWGMSSMNeutral/mssm_xs_tools.py_v2.1
wget -P CombineHarvester/MSSMvsSMRun2Legacy/src -O mssm_xs_tools.cc https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWGMSSMNeutral/mssm_xs_tools.C_v2.1
wget -P CombineHarvester/MSSMvsSMRun2Legacy/interface -O mssm_xs_tools.h https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWGMSSMNeutral/mssm_xs_tools.h_v2.1

# Compile everything
scramv1 b clean; scramv1 b -j $NUM_CORES

# Download root files for latest MSSM benchmark scenarios
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWGMSSMNeutral/mh125_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWGMSSMNeutral/mh125_ls_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWGMSSMNeutral/mh125_lc_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWGMSSMNeutral/mh125_align_13.root
