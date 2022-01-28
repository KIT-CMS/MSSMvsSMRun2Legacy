#!/bin/bash
set -e

# Download old and new model root files from the twiki page and zenodo, respectively.
for model in hMSSM mHH125 mh1125_CPV mh125EFT mh125EFT_lc mh125 mh125_align mh125_lc mh125_ls mh125_muneg_1 mh125_muneg_2 mh125_muneg_3;
do
    wget -O CombineHarvester/MSSMvsSMRun2Legacy/data/${model}_13.root https://zenodo.org/record/5730271/files/${model}_13.root?download=1;
    if [[ "${model}" != "hMSSM" ]]; then
    wget -O CombineHarvester/MSSMvsSMRun2Legacy/data/${model}_13_old.root https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHWGMSSMNeutral/${model}_13.root
    fi
done;
