#!/bin/bash

export SCRAM_ARCH=slc7_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# pushd CMSSW_10_2_25/src
eval `scramv1 runtime -sh`
# popd
