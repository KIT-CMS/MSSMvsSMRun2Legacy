#!/bin/bash
ulimit -s unlimited

INPUT=$1
MODE=$2
CHANNEL=$3

if [[ "$CHANNEL" == "mt" ]]; then
    var="mt_1_puppi"
elif [[ "$CHANNEL" == "em" ]]; then
    var="pzeta"
else
    echo "[ERROR] Given channel not implemented..."
    exit 1
fi

datacarddir=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/categorisation-plots-output/datacards

case "$MODE" in
    "datacards")
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
        [[ ! -d ${datacarddir}/combined ]] && mkdir -p ${datacarddir}/combined
        rsync -av --progress ${datacarddir}/201?/htt_${CHANNEL}_30*/* ${datacarddir}/combined
        ;;

    "ws")
        combineTool.py -M T2W \
            -o ws.root \
            -i ${datacarddir}/{2016,2017,2018}/htt_*/ \
            --parallel 4
        combineTool.py -M T2W \
            -o ws.root \
            -i ${datacarddir}/combined/
        ;;

    "prefit-shapes")
        prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/201*/htt_*/combined.txt.cmb" \
                                          --workspace_name ws.root \
                                          --output_name control-datacard-shapes-prefit.root \
                                          --parallel 8 |& tee -a ${datacarddir}/prefit-shapes-creation.log

        hadd -f ${datacarddir}/combined/control-datacard-shapes-prefit.root ${datacarddir}/201?/htt_*/control-datacard-shapes-prefit.root |& tee -a ${datacarddir}/prefit-shapes-creation.log
        # PostFitShapesFromWorkspace \
        #     -w ${datacarddir}/combined/ws.root \
        #     -d ${datacarddir}/combined/combined.txt.cmb \
        #     -o ${datacarddir}/combined/control-datacard-shapes-prefit.root \
        #     -m 400
        ;;

    "ml-fit")
        combineTool.py -M FitDiagnostics \
            -d ${datacarddir}/combined/ws.root \
            -m 400 \
            --setParameters r=0 --setParameterRange r=-2,2 \
            --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
            --robustHesse 1 \
            -n .combined\
            --there \
            -v 1 |& tee ${datacarddir}/combined/fit_diagnositcs.log
        ;;

    "postfit-shapes")
        PostFitShapesFromWorkspace \
            -w ${datacarddir}/combined/ws.root \
            -d ${datacarddir}/combined/combined.txt.cmb \
            -o ${datacarddir}/combined/control-datacard-shapes-postfit-b.root \
            -f ${datacarddir}/combined/fitDiagnostics.combined.root:fit_b --postfit --sampling --skip-prefit \
            -m 400
        ;;

    "plots")
        source utils/setup_python.sh
        categories="None"
        for FILE in "${datacarddir}/combined/control-datacard-shapes-prefit.root" "${datacarddir}/combined/control-datacard-shapes-postfit-b.root"
        do
            for OPTION in "" "--png"
            do
                for era in 2016 2017 2018 combined; do
                    ./plotting/plot_shapes_gof_categories.py -i $FILE -c $CHANNEL -e $era $OPTION \
                        --categories ${categories} --fake-factor --embedding \
                        --gof-variable ${var} -o ${datacarddir} --linear
                done
            done
        done
        ;;

    *)
        echo "[ERROR] Given mode $MODE is not defined. Aborting..."
        exit 1
        ;;
esac
