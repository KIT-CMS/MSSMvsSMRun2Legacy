#!/bin/bash

source utils/setup_cvmfs_sft.sh
source utils/setup_python.sh

ERA=$1
INPUT=$2
OUTPUT=$3  # Output directory the plots are written to.
IFS="," read -a CHANNELS <<< $4
MASS=$5
XSEC_GGH=$6
XSEC_BBH=$7

[[ -z $MASS ]] && MASS=1200
[[ -z $XSEC_GGH ]] && XSEC_GGH=0.05
[[ -z $XSEC_BBH && ! -z $6 ]] && XSEC_BBH=0.05
[[ -z $XSEC_BBH ]] && XSEC_BBH=0.05


if [[ ! -d "$OUTPUT" ]]
then
    mkdir $OUTPUT
fi

for FILE in $INPUT
do
    for OPTION in "" "--png"
    do
        ./plotting/plot_shapes_mssm.py \
                                       -i $FILE \
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
                                       --x-sec-ggh $XSEC_GGH \
                                       --x-sec-bbh $XSEC_BBH \
    done
done
