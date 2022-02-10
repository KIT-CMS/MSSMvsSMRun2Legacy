#!/bin/bash

ERA=$1
CHANNEL=$2
CATEGORY=$3
DIR=$4
MODE=$5
# Set default value for running mode and check given values for running mode.
[[ -z $5 ]] && MODE="submit"
if [[ "$MODE" != "local" && "$MODE" != "submit" ]]; then
    echo "[ERROR] Given mode $MODE not known. Aborting..."
    exit 1
fi

if [[ "${ERA}" == "combined" ]]; then
    ID=${ERA}
elif [[ "${CHANNEL}" == "cmb" ]]; then
    ID=${ERA}
elif [[ "${CATEGORY}" == "all" ]]; then
    ID=${ERA}-${CHANNEL}
else
    ID=${ERA}-${CHANNEL}-${CATEGORY}
fi

# Switch to working directory, necessary for jobs to run at the right place
# and with the correct environment.
pushd $DIR

# Setup the correct environment for the jobs to run on the node.
source utils/setup_cmssw.sh
# Increase stack size for large workspaces.
ulimit -s unlimited
# Get the code to calculate the correct masks for the given options.
source ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/gof/build_masks.sh

DATACARD=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/analysis/cmb_ind_unblinding/datacards_bsm-model-indep/${ERA}/${CHANNEL}/ws-gof.root
if [[ ! -f ${DATACARD} ]]; then
    echo "Workspace ${DATACARD} does not exist. Please produce it first."
    exit 1
fi

# Change to datacard dir to run fits.
pushd $(dirname $DATACARD)

MASS=160
NUM_TOYS=50 # multiply x10

if [[ ! -d $(dirname ${DATACARD})/gof/${ID}/${ERA}_plots ]]
then
    mkdir -p $(dirname ${DATACARD})/gof/${ID}/${ERA}_plots
fi

MASKS=$(build_masks $ERA $CHANNEL $CATEGORY mod-indep)
MASKS_EVAL=$(build_masks_evaluation $ERA $CHANNEL $CATEGORY mod-indep)
MASK_ARG="--setParametersForFit $MASKS --setParametersForEval ${MASKS_EVAL}"

case "$ERA" in
    "2016")
        TITLE="36.3 fb^{-1} (2016, 13 TeV)"
        ;;
    "2017")
        TITLE="41.5 fb^{-1} (2017, 13 TeV)"
        ;;
    "2018")
        TITLE="59.7 fb^{-1} (2018, 13 TeV)"
        ;;
    "combined")
        TITLE="138 fb^{-1} (13 TeV)"
        ;;
    *)
        echo "[ERROR] Given era $ERA is not defined. Aborting..."
        exit 1
        ;;
esac

for ALGO in "saturated" # "KS" "AD"
do
    # Get test statistic value
    if [[ "$ALGO" == "saturated" ]]
    then
        combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -v 1 --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
    else
        combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -v 1 --setParameters $MASKS --fixedSignalStrength=0 --plots --expectSignal=0 $MASK_ARG
    fi

    # Throw toys
    TOYSOPT=""
    [[ "$ALGO" == "saturated" ]] && TOYSOPT="--toysFreq"

    case "$MODE" in

        local)
            combineTool.py -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1230:1239:1 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --parallel 10 --expectSignal=0 $MASK_ARG
            ;;

        submit)
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1230 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1231 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1232 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1233 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1234 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1235 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1236 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1237 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1238 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            combine -M GoodnessOfFit -n Test.${ID}.${ALGO} --algo=$ALGO -m $MASS -d $DATACARD -s 1239 -t $NUM_TOYS $TOYSOPT --setParameters $MASKS --fixedSignalStrength=0 --expectSignal=0 $MASK_ARG
            ;;
    esac

    # Collect results
    combineTool.py -M CollectGoodnessOfFit --input \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1230.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1231.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1232.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1233.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1234.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1235.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1236.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1237.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1238.root \
        higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.1239.root \
        --output $(dirname ${DATACARD})/gof/${ID}/gof_${ALGO}.json

    if [[ "$ALGO" == "saturated" ]]
    then
        mv $(dirname ${DATACARD})/gof/${ID}/gof_${ALGO}.json $(dirname ${DATACARD})/gof/${ID}/gof.json
    fi

    # Plot
    if [[ "$ALGO" != "saturated" ]]
    then
        plotGof.py --statistic $ALGO --mass $MASS.0 --output gof_${ALGO} $(dirname ${DATACARD})/gof/${ID}/gof_${ALGO}.json --title-right="$TITLE" --title-left="${CHANNEL}, ${CATEGORY}"
        mv htt_*gof_${ALGO}.p{df,ng} $(dirname ${DATACARD})/gof/${ID}/
        ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/gof/plot_gof_metrics.py -e $ERA -g $ALGO -o $(dirname ${DATACARD})/gof/${ID}/${ERA}_plots -i higgsCombineTest.${ID}.${ALGO}.GoodnessOfFit.mH$MASS.root
    else
        plotGof.py --statistic $ALGO --mass $MASS.0 --output gof $(dirname ${DATACARD})/gof/${ID}/gof.json --title-right="$TITLE" --title-left="${CHANNEL}, ${CATEGORY}"
        mv gof.p{df,ng} $(dirname ${DATACARD})/gof/${ID}/
    fi
done
