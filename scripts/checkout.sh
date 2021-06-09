#!/bin/bash

### Setup of CMSSW release
NUM_CORES=10
CMSSW=CMSSW_10_2_25

export SCRAM_ARCH=slc7_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project $CMSSW; pushd $CMSSW/src
eval `scramv1 runtime -sh`

# combine on 102X slc7
#git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
#cd HiggsAnalysis/CombinedLimit
#git fetch origin
#git checkout v8.0.1
#cd -
git clone git@github.com:KIT-CMS/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

# CombineHarvester (current master for 102X)
#git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
git clone git@github.com:KIT-CMS/CombineHarvester.git CombineHarvester # fixes & extensions for morphing (22.08.2019)

# MSSM vs SM analysis specific code
git clone git@github.com:KIT-CMS/MSSMvsSMRun2Legacy.git CombineHarvester/MSSMvsSMRun2Legacy

# SM analysis specific code
git clone git@github.com:KIT-CMS/SMRun2Legacy.git CombineHarvester/SMRun2Legacy

# grid control
git clone git://github.com/harrypuuter/grid-control.git

# Install LHCHWGMSSMNeutral packages
wget -O CombineHarvester/MSSMvsSMRun2Legacy/python/mssm_xs_tools.py https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mssm_xs_tools.py_v2.4
sed -i "s!./mssm_xs_tools_C.so!${CMSSW_BASE}/lib/${SCRAM_ARCH}/libCombineHarvesterMSSMvsSMRun2Legacy.so!g" CombineHarvester/MSSMvsSMRun2Legacy/python/mssm_xs_tools.py
sed -i "s!root_files/mh125_13.root!${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/mh125_13.root!g" CombineHarvester/MSSMvsSMRun2Legacy/python/mssm_xs_tools.py
wget -O CombineHarvester/MSSMvsSMRun2Legacy/src/mssm_xs_tools.cc https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mssm_xs_tools.C_v2.4
sed -i 's!mssm_xs_tools.h!CombineHarvester/MSSMvsSMRun2Legacy/interface/mssm_xs_tools.h!g' CombineHarvester/MSSMvsSMRun2Legacy/src/mssm_xs_tools.cc
wget -O CombineHarvester/MSSMvsSMRun2Legacy/interface/mssm_xs_tools.h https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mssm_xs_tools.h_v2.4

# Compile everything
scramv1 b clean; scramv1 b -j $NUM_CORES

# Download root files for latest MSSM benchmark scenarios
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125_ls_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125_lc_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125_align_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mHH125_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh1125_CPV_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125_muneg_1_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125_muneg_2_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125_muneg_3_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125EFT_13.root
wget -P CombineHarvester/MSSMvsSMRun2Legacy/data https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/mh125EFT_lc_13.root

# Download ggH NLO reweighting inputs
wget https://github.com/danielwinterbottom/ggh-mssm/raw/master/workspace/higgs_pt_v2.root -O CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root
