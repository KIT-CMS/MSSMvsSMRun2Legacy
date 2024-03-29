#!/bin/bash
ulimit -s unlimited

INPUT=$1
MODE=$2
CHANNEL=$3

if [[ "$CHANNEL" == mt ]]; then
    var="mt_1_puppi"
elif [[ "$CHANNEL" == em ]]; then
    var="pzeta"
elif [[ "$CHANNEL" == cmb ]]; then
    var=""
    if [[ "$MODE" == "datacards" ]]; then
        echo "[ERROR] Datacard step not defined for channel $CHANNEL"
        exit 1
    fi
else
    echo "[ERROR] Given channel not implemented..."
    exit 1
fi

datacarddir=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/categorisation-plots-output/datacards

case "$MODE" in
    "datacards")
        if [[ ${CHANNEL} != mt && ${CHANNEL} != em ]]; then
            echo "[ERROR] Channel for datacard step must be either em or mt."
            exit 1
        fi
        [[ ! -d  ${datacarddir} ]] && mkdir -p ${datacarddir}
        for era in 2016 2017 2018; do
            MorphingCatVariables \
                --base-path ${INPUT} \
                --output_folder=${datacarddir} \
                --era=${era} \
                --category=all  --variable=${var} \
                --mode=categorisation-plots \
                --channel $CHANNEL \
                --use_mc=0 \
                --verbose=1
        done
        [[ ! -d ${datacarddir}/combined/cmb ]] && mkdir -p ${datacarddir}/combined/cmb
        rsync -av --progress ${datacarddir}/201?/htt_${CHANNEL}_30*/* ${datacarddir}/combined/cmb
        [[ ! -d ${datacarddir}/combined/${CHANNEL} ]] && mkdir -p ${datacarddir}/combined/${CHANNEL}
        rsync -av --progress ${datacarddir}/201?/htt_${CHANNEL}_30*/* ${datacarddir}/combined/${CHANNEL}
        [[ ! -d ${datacarddir}/combined/htt_${CHANNEL}_301 ]] && mkdir -p ${datacarddir}/combined/htt_${CHANNEL}_301
        rsync -av --progress ${datacarddir}/201?/htt_${CHANNEL}_301*/* ${datacarddir}/combined/htt_${CHANNEL}_301
        [[ ! -d ${datacarddir}/combined/htt_${CHANNEL}_302 ]] && mkdir -p ${datacarddir}/combined/htt_${CHANNEL}_302
        rsync -av --progress ${datacarddir}/201?/htt_${CHANNEL}_302*/* ${datacarddir}/combined/htt_${CHANNEL}_302

        # Create workspaces for categories directly when creating datacards
        combineTool.py -M T2W \
            -o ws.root \
            -i ${datacarddir}/{2016,2017,2018,combined}/htt_*/ \
            --parallel 4
        ;;

    "ws")
        combineTool.py -M T2W \
            -o ws.root \
            -i ${datacarddir}/combined/${CHANNEL}/
        ;;

    "prefit-shapes")
        prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/201*/htt_*/combined.txt.cmb" \
                                          --workspace_name ws.root \
                                          --output_name control-datacard-shapes-prefit.root \
                                          --parallel 8 |& tee ${datacarddir}/prefit-shapes-creation.log

        hadd -f ${datacarddir}/combined/cmb/control-datacard-shapes-prefit.root ${datacarddir}/201?/htt_*/control-datacard-shapes-prefit.root |& tee -a ${datacarddir}/prefit-shapes-creation.log
        for ch in mt em; do
            hadd -f ${datacarddir}/combined/${ch}/control-datacard-shapes-prefit.root ${datacarddir}/201?/htt_${ch}*/control-datacard-shapes-prefit.root |& tee -a ${datacarddir}/prefit-shapes-creation.log
        done
        ;;

    "ml-fit")
        combineTool.py -M FitDiagnostics \
            -d ${datacarddir}/combined/${CHANNEL}/ws.root \
            -m 400 \
            --setParameters r=0 --setParameterRange r=-2,2 \
            --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
            --robustHesse 1 \
            -n .combined.${CHANNEL} \
            --there \
            -v 1 |& tee ${datacarddir}/combined/${CHANNEL}/fit_diagnositcs.log
        ;;

    "postfit-shapes")
        ch_flag=""
        [[ "$CHANNEL" != "cmb" ]] && ch_flag=${CHANNEL}
        PostFitShapesFromWorkspace \
            -d ${datacarddir}/combined/${CHANNEL}/combined.txt.cmb \
            -w ${datacarddir}/combined/${CHANNEL}/ws.root \
            -o ${datacarddir}/combined/${CHANNEL}/control-datacard-shapes-postfit-b-${CHANNEL}_fit.root \
            -m 400 \
            -f ${datacarddir}/combined/${fit_ch}/fitDiagnostics.combined.${fit_ch}.root:fit_b \
            --postfit --sampling --skip-prefit --total-shapes
        for cat in 301 302; do
            PostFitShapesFromWorkspace \
                -d ${datacarddir}/combined/htt_${CHANNEL}_${cat}/combined.txt.cmb \
                -w ${datacarddir}/combined/htt_${CHANNEL}_${cat}/ws.root \
                -o ${datacarddir}/combined/htt_${CHANNEL}_${cat}/control-datacard-shapes-postfit-b-${CHANNEL}_fit.root \
                -m 400 \
                -f ${datacarddir}/combined/${CHANNEL}/fitDiagnostics.combined.${CHANNEL}.root:fit_b \
                --postfit --sampling --skip-prefit --total-shapes
        done
        ;;

    "plots")
        source utils/setup_python.sh
        [[ ! -d ${datacarddir}/plots/${CHANNEL}-fit ]] && mkdir -p ${datacarddir}/plots/${CHANNEL}-fit
        categories="None"
        for era in 2016 2017 2018; do
            for FILE in "${datacarddir}/combined/${CHANNEL}/control-datacard-shapes-prefit.root" "${datacarddir}/combined/${CHANNEL}/control-datacard-shapes-postfit-b-${CHANNEL}_fit.root"
            do
                for OPTION in "" "--png"
                do
                    ./plotting/plot_shapes_gof_categories.py -i $FILE -c $CHANNEL -e $era $OPTION \
                        --categories ${categories} --fake-factor --embedding \
                        --gof-variable ${var} -o ${datacarddir}/plots/${CHANNEL}-fit --linear
                done
            done
        done

        # Run combined separately as special treatment for postfit is needed
        for FILE in "${datacarddir}/combined/${CHANNEL}/control-datacard-shapes-prefit.root" "${datacarddir}/combined/${CHANNEL}/control-datacard-shapes-postfit-b-${CHANNEL}_fit.root" \
                    "${datacarddir}/combined/htt_${CHANNEL}_301/control-datacard-shapes-postfit-b-${CHANNEL}_fit.root" "${datacarddir}/combined/htt_${CHANNEL}_302/control-datacard-shapes-postfit-b-${CHANNEL}_fit.root"
        do
            for OPTION in "" "--png"
            do
                ./plotting/plot_shapes_gof_categories.py -i $FILE -c $CHANNEL -e combined $OPTION \
                    --categories ${categories} --fake-factor --embedding \
                    --gof-variable ${var} -o ${datacarddir}/plots/${CHANNEL}-fit --linear
            done
        done
        ;;

    *)
        echo "[ERROR] Given mode $MODE is not defined. Aborting..."
        exit 1
        ;;
esac
