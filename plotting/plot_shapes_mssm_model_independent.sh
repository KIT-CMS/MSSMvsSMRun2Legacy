#!/bin/bash

source utils/setup_cvmfs_sft.sh
source utils/setup_python.sh

ERA=$1
INPUT=$2
OUTPUT=$3  # Output directory the plots are written to.
IFS="," read -a CHANNELS <<< $4
MASS=$5
XSEC=$6

[[ -z $MASS ]] && MASS=1200
[[ -z $XSEC ]] && XSEC=0.05


if [[ ! -d "$OUTPUT" ]]
then
    mkdir $OUTPUT
fi

for FILE in $INPUT
do
    for OPTION in "" "--png"
    do
        ./plotting/plot_shapes_mssm.py -i $FILE \
                                       -c ${CHANNELS[@]} \
                                       -e $ERA \
                                       $OPTION \
                                       --fake-factor \
                                       --embedding \
                                       --normalize-by-bin-width \
                                       -o $OUTPUT \
                                       --model-independent \
                                       --blinded \
                                       --mass $MASS \
                                       --x-sec $XSEC
    done
done
