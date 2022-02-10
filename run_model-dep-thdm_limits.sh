#!/bin/bash

ulimit -s unlimited

TAG=$1
MODE=$2
MODEL=$3
ANALYSISTYPE=$4

analysis="bsm-model-dep-additional"
sm_like_hists="sm125"
replace_with_sm125=1
if [[ "$ANALYSISTYPE" == "classic" ]]; then
    categorization="classic"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_classic"
    fi
else
    categorization="with-sm-ml"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_with_ml"
    fi
fi

case "$MODEL" in
    "THDM_BP1_TYPE1")
        wsoutput="ws_thdm_bp1_type1.root"
        modelfile="BP1_Type1.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='m_{H} [GeV]'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    "THDM_BP1_TYPE2")
        wsoutput="ws_thdm_bp1_type2.root"
        modelfile="BP1_Type2.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='m_{H} [GeV]'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    "THDM_FixedMass_TYPE2")
        wsoutput="ws_thdm_fixedmass_type2.root"
        modelfile="FixedMass_Type2.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='cos(#beta-#alpha)'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    *)
        echo "[ERROR] Given model $MODEL not known..."
        exit 1
esac

defaultdir="analysis_v14/$TAG"
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
defaultdir=$(readlink -f ${defaultdir})
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/limits_${MODEL}/condor ]] && mkdir -p ${defaultdir}/limits_${MODEL}/condor

datacarddir=${defaultdir}/datacards_${analysis}
taskname="${analysis}_${TAG}_${MODEL}_1"
taskname2="${analysis}_${TAG}_${MODEL}_2"

if [[ $MODE == "initial" ]]; then
    ############
    # morphing
    ############
    if [[ $ANALYSISTYPE == "classic" ]]; then
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment "hSM-in-bg" \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=0 --manual_rebin=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
    else
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment "hSM-in-bg" \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=0 --manual_rebin=1 --split_sm_signal_cat=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_neuralnet_categories.txt \
            --variable nnscore \
            --sm \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment "hSM-in-bg" \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=0 --manual_rebin=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_new_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
    fi

    ############
    # combining outputs
    ############
    mkdir -p ${datacarddir}/combined/cmb/

    rsync -av --progress ${datacarddir}/201?/htt_*/* ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
    for era in 2016 2017 2018;
    do
        mkdir -p ${datacarddir}/${era}/cmb/
        rsync -av --progress ${datacarddir}/${era}/htt_*/* ${datacarddir}/${era}/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards_${era}.txt
        for channel in "et" "mt" "tt";
        do
            mkdir -p ${datacarddir}/${era}/${channel}/
            rsync -av --progress ${datacarddir}/${era}/htt_${channel}*/* ${datacarddir}/${era}/${channel}/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards_${era}_${channel}.txt
        done
    done

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############
    combineTool.py -M T2W -o ${wsoutput} \
    -P CombineHarvester.MSSMvsSMRun2Legacy.THDMvsSM:THDMvsSM \
    --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
    --PO modelFile=${modelfile} \
    --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
    -i ${datacarddir}/2016/mt/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt

elif [[ "$MODE" == "setup" ]]; then
    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/thdm_asymptotic_grid_${MODEL}.json \
    -d ${datacarddir}/2016/mt/${wsoutput} \
    --job-mode 'condor' \
    --task-name $taskname \
    --dry-run \
    --redefineSignalPOI r \
    --setParameterRanges r=0,1 \
    --setParameters r=1,x=1 \
    --freezeParameters x -v 1 \
    -t -1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --merge 5 \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}.txt

elif [[ $MODE == "submit" ]]; then
    ############
    # job submission
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    condor_submit condor_${taskname}.sub

elif [[ $MODE == "submit-local" ]]; then
    ############
    # job submission
    ############
    cp scripts/run_limits_locally.py ${defaultdir}/limits_${MODEL}/condor
    cd ${defaultdir}/limits_${MODEL}/condor
    python run_limits_locally.py --cores 20 --taskname condor_${taskname}.sh

elif [[ $MODE == "collect" ]]; then
    ############
    # job collection
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/thdm_asymptotic_grid_${MODEL}.json \
    -d ${datacarddir}/2016/mt/${wsoutput} \
    --job-mode 'condor' \
    --task-name $taskname2 \
    --dry-run \
    --redefineSignalPOI x \
    --setParameterRanges x=0,1 \
    --setParameters r=1 \
    --freezeParameters r -v 1 \
    -t -1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/collect_jobs_${MODEL}.txt

    # # condor_submit condor_${taskname2}.sub
    cp asymptotic_grid.root ..
    cd ${defaultdir}/limits_${MODEL}/

    ############
    # limit plot
    ############
    if [[ $ANALYSISTYPE == "classic" ]]; then
        title="Classic categorisation 36 fb^{-1} (2016, 13 TeV)"
    else
        title="138 fb^{-1} (13 TeV)"
    fi
    plotLimitGrid.py asymptotic_grid.root \
    --scenario-label="${scenario_label}" \
    --output ${TAG}_${MODEL} \
    --title-right="${title}" \
    --cms-sub="Preliminary" \
    --contours="exp-2,exp-1,exp0,exp+1,exp+2" \
    --thdm-validity --logy \
    --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${modelfile} \
    --x-title "${x_title}" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt
    # --y-range ${y_min},${y_max} \
    # --mass_histogram ${sm_like_mass} \
    # --mass_histogram_title ${mass_histogram_title} \
fi
