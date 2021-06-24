#!/bin/bash

# Ensure stack size is large enough
ulimit -s unlimited

TAG=$1
MODE=$2
MASS=$3
PROCESS=$4
ERA=$5
CHANNEL=$6

# Sanity checks on input
[[ ! -z $1 && ! -z $2 && ! -z $3 && ! -z $4 ]] || ( echo "[ERROR] Number of given parameters is too small."; exit 1 )
[[ -z $5 ]] && echo "[INFO] Era not set. Will use combined dataset." && ERA="combined"
[[ -z $6 ]] && echo "[INFO] Channel not set. Will use combination of channels." && CHANNEL="cmb"

if [[ "${TAG}" == "auto" ]]; then
    TAG=${ERA}_${CHANNEL}
fi

defaultdir="${CMSSW_BASE}/src/fits_powheg/${TAG}"
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
defaultdir=$(readlink -f ${defaultdir})
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/impacts_${PROCESS}_mH${MASS}/condor ]] && mkdir -p ${defaultdir}/impacts_${PROCESS}_mH${MASS}/condor
if [[ "$PROCESS" == "r_ggH" ]]; then
    signal_process=$PROCESS
    freeze_process=r_bbH
elif [[ "$PROCESS" == "r_bbH" ]]; then
    signal_process=$PROCESS
    freeze_process=r_ggH
else
    echo "[ERROR] Process not known aborting."
    exit 1
fi


datacarddir=$(dirname ${defaultdir})/output_mssm_classic
taskname="impacts_${TAG}_${PROCESS}_mH${MASS}"

case $MODE in

    prepare-submit)
        pushd ${defaultdir}/impacts_${PROCESS}_mH${MASS}/condor
        combineTool.py -M Impacts -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
                       --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
                       --doInitialFit --robustFit 1 \
                       -t -1 -m $MASS \
                       --setParameters r_ggH=0,r_bbH=0 --setParameterRanges r_ggH=-1.0,1.0:r_bbH=-1.0,1.0 \
                       --redefineSignalPOIs ${signal_process} --freezeParameters ${freeze_process} -v 1 # --stepSize 0.01
                       # --rAbsAcc 0 --rRelAcc 0.0005 \

        combineTool.py -M Impacts -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
                       --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
                       --robustFit 1 --doFits \
                       -t -1 -m $MASS \
                       --setParameters r_ggH=0,r_bbH=0 --setParameterRanges r_ggH=-1.0,1.0:r_bbH=-1.0,1.0 \
                       --redefineSignalPOIs ${signal_process} --freezeParameters ${freeze_process} \
                       --job-mode condor --task-name ${taskname} --dry-run --merge 5
                       # --rAbsAcc 0 --rRelAcc 0.0005 \
        popd
        ;;

    submit)
        ############
        # job submission
        ############
        pushd ${defaultdir}/impacts_${PROCESS}_mH${MASS}/condor
        condor_submit condor_${taskname}.sub
        popd
        ;;

    submit-gc)
        ############
        # job submission
        ############
        python scripts/build_gc_job.py \
            --combine-script ${defaultdir}/impacts_${PROCESS}_mH${MASS}/condor/condor_${taskname}.sh \
            --workspace ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
            --workdir /work/mburkart/workdirs/combine/${taskname} \
            --tag ${taskname} \
            --se-path /storage/gridka-nrg/mburkart/gc_storage/combine/${TAG}/${taskname}

        ${CMSSW_BASE}/src/grid-control/go.py /work/mburkart/workdirs/combine/${taskname}/${taskname}.conf -Gc -m 3
        ;;

    copy-results-gc)
        ############
        # job submission
        ############
        rsync -avhP /storage/gridka-nrg/mburkart/gc_storage/combine/${TAG}/${taskname}/output/ ${defaultdir}/impacts_${PROCESS}_mH${MASS}/condor
        ;;

    collect)
        ############
        # job collection
        ############
        pushd ${defaultdir}/impacts_${PROCESS}_mH${MASS}/condor
        combineTool.py -M Impacts -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root -m $MASS -o ${ERA}_${CHANNEL}_${signal_process}_${MASS}_impacts.json --redefineSignalPOIs ${signal_process} --exclude ${freeze_process}

        plotImpacts.py -i ${ERA}_${CHANNEL}_${signal_process}_${MASS}_impacts.json -o ${ERA}_${CHANNEL}_${signal_process}_${MASS}_impacts --transparent --translate ${CMSSW_BASE}/src/translate.json # _nobbb --no-bbb
        popd
        ;;

    *)
        echo "[ERROR] Given mode not known. Skipping."
        exit 1
        ;;
esac
