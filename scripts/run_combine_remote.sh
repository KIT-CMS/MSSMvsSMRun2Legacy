#!/bin/bash
ulimit -s unlimited
set -e
export XRD_LOGLEVEL="Info"
export SCRAM_ARCH=${SCRAM_ARCH}
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $GC_GLITE_LOCATION
echo " --------------"
echo "Job Parameter Set: "
echo "Running on " hostname
echo "TARBALL_PATH: ${SCRAM_ARCH}"
echo "TARBALL_PATH: ${TARBALL_PATH}"
echo "OUTPUT_PATH: ${OUTPUT_PATH}"
echo "CMSSW_VERSION: ${CMSSW_VERSION}"
echo "GC_GLITE_LOCATION: ${GC_GLITE_LOCATION}"
echo "JOB_ID: ${JOB_ID}"
echo "--------------"
echo "setting up cmssw"
gfal-copy ${TARBALL_PATH} .
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project ${CMSSW_VERSION}
tar -xf cmssw.tar.gz
rm cmssw.tar.gz
cd ${CMSSW_VERSION}/src
scram b ProjectRename
eval $(scram runtime -sh)
cd -

echo " --------------"
echo " Starting Combine"
./combine_fit.sh ${JOB_ID}
rm ws.root


echo " --------------"
echo " Finished Fit !"
echo xrdcp -pf *.root ${OUTPUT_PATH}
xrdcp -pf *.root ${OUTPUT_PATH}

echo " --------------"
echo " Finished Fit !"