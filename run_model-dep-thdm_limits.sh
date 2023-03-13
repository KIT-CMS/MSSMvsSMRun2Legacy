#!/bin/bash

ulimit -s unlimited

TAG=$1
MODE=$2
MODEL=$3
ANALYSISTYPE=$4

GRIDUSER=$(whoami)
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
    "THDM_BP1_Type1")
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
    "THDM_BP1_Type1_mod")
        wsoutput="ws_thdm_bp1_type1_mod.root"
        modelfile="BP1_Type1_mod.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='m_{H} [GeV]'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    "THDM_BP1_Type2")
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
    "THDM_FixedMass_Type2")
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
    "THDM_PASComp_Type1")
        wsoutput="ws_thdm_pascomp_type1.root"
        modelfile="PASComp_Type1.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='m_{H} [GeV]'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    "THDM_PASComp_Type2")
        wsoutput="ws_thdm_pascomp_type2.root"
        modelfile="PASComp_Type2.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='m_{H} [GeV]'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    "THDM_HWWLike_Type1")
        wsoutput="ws_thdm_hwwlike_type1.root"
        modelfile="HWWLike_Type1.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='m_{H} [GeV]'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    "THDM_HWWLike_Type2")
        wsoutput="ws_thdm_hwwlike_type2.root"
        modelfile="HWWLike_Type2.root"
        # scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
        sub_analysis="sm-like-light"
        # sm_like_mass="m_h"
        x_title='m_{H} [GeV]'
        # mass_histogram_title="m_{h}"
        # y_min=1.0
        # y_max=60.0
        ;;
    *)
        echo -e "\033[0;31m[ERROR]\033[0m Given model $MODEL not known..."
        exit 1
esac

defaultdir="analysis_2022_02_28/$TAG"
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
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
    elif [[ $ANALYSISTYPE == "with-sm-ml" ]]; then
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment "hSM-in-bg" \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1" \
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
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_new_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
    else
        echo -e "\033[0;31m[ERROR]\033[0m Given analysis type `${ANALYSISTYPE}` not known. Please try again..."
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
    # Check if the expected number of datacards has been written
    if [[ $ANALYSISTYPE == "classic" ]]; then
        EXPECTED=$(((4+4+2+7)*3))
    else
        EXPECTED=$(((15+15+11+18)*3))
    fi
    if [[ $(ls ${datacarddir}/combined/cmb/*.txt | wc -l) != $EXPECTED ]]; then
        echo -e "\033[0;31m[ERROR]\033[0m Not all datacards have been created or written. Please check the logs..."
        echo "Expected ${EXPECTED} datacards written but found only $(ls ${datacarddir}/combined/cmb/*.txt | wc -l) in the combined directory."
    fi

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############
    combineTool.py -M T2W -o ${wsoutput} \
    -P CombineHarvester.MSSMvsSMRun2Legacy.THDMvsSM:THDMvsSM \
    --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
    --PO modelFile=${modelfile} \
    --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
    --PO replace-with-sm125=${replace_with_sm125} \
    --PO sm-predictions=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_predictions_13TeV.json \
    -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt
    # Not needed as only relevant for full limits including h signal hypothesis
    # --PO hSM_treatment="hSM-in-bg" \

elif [[ "$MODE" == "setup" ]]; then
    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/thdm_asymptotic_grid_${MODEL}.json \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'condor' \
    --task-name $taskname \
    --dry-run \
    --redefineSignalPOI r \
    --setParameterRanges r=0,1 \
    --setParameters r=1,x=1 \
    --freezeParameters x -v 1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --merge 5 \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}.txt
    # -t -1 \

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

elif [[ $MODE == "submit-gc" ]]; then
    ############
    # job submission
    ############
    gcworkdir=${defaultdir}/limits_${MODEL}/gc_condor
    mkdir -p ${gcworkdir}
    python scripts/build_gc_job.py \
        --combine-script ${defaultdir}/limits_${MODEL}/condor/condor_${taskname}.sh \
        --workspace ${datacarddir}/combined/cmb/${wsoutput} \
        --workdir ${gcworkdir} \
        --tag ${taskname} \
        --se-path /storage/gridka-nrg/${GRIDUSER}/gc_storage/combine/${taskname}

    ${CMSSW_BASE}/src/grid-control/go.py ${gcworkdir}/${taskname}.conf -Gc -m 3

elif [[ $MODE == "copy-results-gc" ]]; then
    ############
    # job submission
    ############
    rsync -avhP /storage/gridka-nrg/${GRIDUSER}/gc_storage/combine/${taskname}/output/ ${defaultdir}/limits_${MODEL}/condor

elif [[ $MODE == "collect" ]]; then
    ############
    # job collection
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/thdm_asymptotic_grid_${MODEL}.json \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'condor' \
    --task-name $taskname2 \
    --dry-run \
    --redefineSignalPOI r \
    --setParameterRanges r=0,1 \
    --setParameters x=1,r=1 \
    --freezeParameters x -v 1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/collect_jobs_${MODEL}.txt
    # -t -1 \

    # # condor_submit condor_${taskname2}.sub
    cp asymptotic_grid.root ..
    cd ${defaultdir}/limits_${MODEL}/

    ############
    # limit plot
    ############
    if [[ $ANALYSISTYPE == "classic" ]]; then
        titleleft='--title-left="Classic_categorisation"'
        title="#font[62]{CMS} data 138 fb^{-1} (13 TeV)"
    else
        titleleft=""
        title="#font[62]{CMS} data 138 fb^{-1} (13 TeV)"
    fi
    echo $titleleft
    ${CMSSW_BASE}/src/CombineHarvester/CombineTools/scripts/plotLimitGrid.py asymptotic_grid.root \
    --scenario-label="${scenario_label}" \
    --output ${TAG}_${MODEL}_obs \
    --title-right="${title}" \
    ${titleleft} \
    --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
    --thdm-validity --logy \
    --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${modelfile} \
    --y-title "tan#kern[0.16667]{#beta}" \
    --x-title "${x_title}" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt
    # --y-range ${y_min},${y_max} \
    # --mass_histogram ${sm_like_mass} \
    # --mass_histogram_title ${mass_histogram_title} \
fi
