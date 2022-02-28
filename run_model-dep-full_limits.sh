# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
MODEL=$3
ANALYSISTYPE=$4
HSMTREATMENT=$5
GRIDUSER=$6
CYCLES=$7
OLDFILES=$8
[[ -z $8 ]] && OLDFILES=0

if [[ $ANALYSISTYPE == "classic" ]]; then
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="classic"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_classic"
    fi
elif [[ $ANALYSISTYPE == "classic_lowmass" ]]; then
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="lowmass"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_classic_lowmass"
    fi
elif [[ $ANALYSISTYPE == "with-sm-ml" ]]; then
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="with-sm-ml"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_with_ml"
    fi
elif [[ $ANALYSISTYPE == "sm-ml-only" ]]; then
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="sm-ml-only"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_sm_ml_only"
    fi
else
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="with-sm-ml"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_with_ml"
    fi
fi
# Szenarios from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGMSSMNeutral?redirectedfrom=LHCPhysics.LHCHXSWGMSSMNeutral#Baseline_scenarios
### MSSM scenarios #####
scale_qqh_by_hand=0
if [[ $MODEL == "mh125" ]]; then
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    scenario_label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_lc" ]]; then
    wsoutput="ws_mh125_lc.root"
    modelfile="13,Run2017,mh125_lc_13.root"
    scenario_label="M_{h}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_ls" ]]; then
    wsoutput="ws_mh125_ls.root"
    modelfile="13,Run2017,mh125_ls_13.root"
    scenario_label="M_{h}^{125}(#tilde{#tau}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_align" ]]; then
    wsoutput="ws_mh125_align.root"
    modelfile="13,Run2017,mh125_align_13.root"
    scenario_label="M_{h}^{125}(alignment) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=3.0
    y_max=12.0
elif [[ $MODEL == "mHH125" ]]; then
    wsoutput="ws_mHH125.root"
    modelfile="13,Run2017,mHH125_13.root"
    scenario_label="M_{H}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-heavy"
    sm_like_mass="m_H"
    x_title='m_{H^{#plus}} [GeV]'
    mass_histogram_title="m_{H}"
    y_min=5.0
    y_max=6.0
elif [[ $MODEL == "mh1125_CPV" ]]; then
    wsoutput="ws_mh1125_cpv.root"
    modelfile="13,Run2017,mh1125_CPV_13.root"
    scenario_label="M_{h_{1}}^{125}(CPV) scenario (^{}h_{1},^{}h_{2},^{}h_{3}#rightarrow#tau#tau)"
    sub_analysis="cpv"
    sm_like_mass="m_H1"
    x_title='m_{H^{#plus}} [GeV]'
    mass_histogram_title="m_{^{}h_{1}}"
    y_min=1.0
    y_max=20.0
### Negative mu scenarios #####
elif [[ $MODEL == "mh125_muneg_1" ]]; then
    wsoutput="mh125_muneg_1.root"
    modelfile="13,Run2017,mh125_muneg_1_13.root"
    scenario_label="M_{h}^{125 ^{}#mu_{1}#minus} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=4.0
    y_max=56.0
elif [[ $MODEL == "mh125_muneg_2" ]]; then
    wsoutput="mh125_muneg_2.root"
    modelfile="13,Run2017,mh125_muneg_2_13.root"
    scenario_label="M_{h}^{125 ^{}#mu_{2}#minus} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=5.0
    y_max=30.0
elif [[ $MODEL == "mh125_muneg_3" ]]; then
    wsoutput="mh125_muneg_3.root"
    modelfile="13,Run2017,mh125_muneg_3_13.root"
    scenario_label="M_{h}^{125 ^{}#mu_{3}#minus} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=6.0
    y_max=20.0
### EFT scenarios #####
elif [[ $MODEL == "mh125EFT" ]]; then
    wsoutput="mh125EFT.root"
    modelfile="13,Run2017,mh125EFT_13.root"
    scenario_label="M_{h,EFT}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=10.0
elif [[ $MODEL == "mh125EFT_lc" ]]; then
    wsoutput="mh125EFT_lc.root"
    modelfile="13,Run2017,mh125EFT_lc_13.root"
    scenario_label="M_{h,EFT}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=10.0
### hMSSM scenario #####
elif [[ $MODEL == "hMSSM" ]]; then
    wsoutput="hMSSM.root"
    modelfile="13,Run2017,hMSSM_13.root"
    scenario_label="hMSSM scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
else
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    scenario_label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
fi

if [[ $OLDFILES == 1 ]]; then
    wsoutput=${wsoutput/.root/_old.root}
    modelfile=${modelfile/.root/_old.root}
fi

defaultdir="analysis/$TAG"
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
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
    elif [[ $ANALYSISTYPE == "classic_lowmass" ]]; then
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_2d_to_1d.txt \
            --variable m_sv_VS_pt_tt_splitpT \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log_lowmass.txt
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_1d_btag.txt \
            --variable m_sv_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log_btag.txt
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_cr.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log_cr.txt
    elif [[ $ANALYSISTYPE == "with-sm-ml" || $ANALYSISTYPE == "sm-ml-only" ]]; then
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_neuralnet_categories.txt \
            --variable nnscore \
            --sm \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

        if [[ $ANALYSISTYPE == "with-sm-ml" ]]; then
            morph_parallel.py --output ${defaultdir}/datacards \
                --analysis ${analysis} \
                --sub-analysis ${sub_analysis} \
                --hSM-treatment $HSMTREATMENT  \
                --categorization ${categorization} \
                --sm-like-hists ${sm_like_hists} \
                --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
                --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
                --eras 2016,2017,2018 \
                --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_new_categories.txt \
                --variable mt_tot_puppi \
                --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
        fi
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
    EXPECTED=$(((15+15+11+18)*3))
    if [[ $(ls ${datacarddir}/combined/cmb/*.txt | wc -l) != $EXPECTED ]]; then
        echo "[ERROR] Not all datacards have been created or written. Please check the logs..."
        echo "Expected ${EXPECTED} datacards written but found only $(ls ${datacarddir}/combined/cmb/ | wc -l) in the combined directory."
    fi

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############
    if [[ $OLDFILES == 0 ]]; then
        combineTool.py -M T2W -o ${wsoutput} \
        -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM \
        --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
        --PO replace-with-SM125=${replace_with_sm125} \
        --PO hSM-treatment=$HSMTREATMENT \
        --PO modelFile=${modelfile} \
        --PO minTemplateMass=60 \
        --PO maxTemplateMass=3500 \
        --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --PO sm-predictions=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_predictions_13TeV.json \
        --PO qqh-pred-from-scaling=${scale_qqh_by_hand} \
        -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt
    else
        combineTool.py -M T2W -o ${wsoutput} \
        -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM_oldModels:MSSMvsSM_oldModels \
        --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
        --PO replace-with-SM125=${replace_with_sm125} \
        --PO hSM-treatment=$HSMTREATMENT \
        --PO modelFile=${modelfile} \
        --PO minTemplateMass=60 \
        --PO maxTemplateMass=3500 \
        --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --PO sm-predictions=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_predictions_13TeV.json \
        -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt
    fi

elif [[ $MODE == "setup" ]]; then

    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    if [[ $HSMTREATMENT == "hSM-in-bg" ]]; then
        combineTool.py -M AsymptoticGrid \
        ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_${MODEL}.json \
        -d ${datacarddir}/combined/cmb/${wsoutput} \
        --job-mode 'condor' \
        --task-name $taskname \
        --dry-run \
        --redefineSignalPOI r \
        --setParameterRanges r=0,1 \
        --setParameters r=1,x=1 \
        --freezeParameters x -v1 \
        --cminDefaultMinimizerStrategy 0 \
        --X-rtd MINIMIZER_analytic \
        --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}.txt
    elif [[ $HSMTREATMENT == "no-hSM-in-bg" ]]; then
        combineTool.py -M AsymptoticGrid \
        ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_${MODEL}.json \
        -d ${datacarddir}/combined/cmb/${wsoutput} \
        --job-mode 'condor' \
        --task-name $taskname \
        --dry-run \
        --redefineSignalPOI x \
        --setParameterRanges x=0,1 \
        --setParameters r=1 \
        --freezeParameters r -v1 \
        --cminDefaultMinimizerStrategy 0 \
        --X-rtd MINIMIZER_analytic \
        --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}.txt
    fi

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

elif [[ $MODE == "hybrid-lhc" ]]; then

    mkdir -p ${defaultdir}/limits_${MODEL}_hybrid_lhc/condor
    cd ${defaultdir}/limits_${MODEL}_hybrid_lhc/condor

    jsonfile=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_LHC_${MODEL}.json
    if [[ $HSMTREATMENT == "hSM-in-bg" ]]; then
        jsonfile=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_LHC_${MODEL}_newsigmodel.json
    elif [[ $HSMTREATMENT == "no-hSM-in-bg" ]]; then
        jsonfile=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_LHC_${MODEL}.json
    fi

    combineTool.py -M HybridNewGrid \
    ${jsonfile} \
    --cycles $CYCLES \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'crab3' \
    --task-name ${taskname}_hybrid_lhc \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    --X-rtd MINIMIZER_analytic \
    -v1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}_hybrid_lhc.txt
    cd -

elif [[ $MODE == "hybrid-tev" ]]; then

    mkdir -p ${defaultdir}/limits_${MODEL}_hybrid_tev/condor
    cd ${defaultdir}/limits_${MODEL}_hybrid_tev/condor

    combineTool.py -M HybridNewGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_TEV_${MODEL}.json \
    --cycles $CYCLES \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'crab3' \
    --task-name ${taskname}_hybrid_tev \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    --X-rtd MINIMIZER_analytic \
    -v1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}_hybrid_tev.txt
    cd -

elif [[ $MODE == "collect-hybrid-lhc" ]]; then

    cd ${defaultdir}/limits_${MODEL}_hybrid_lhc/condor

    combineTool.py -M HybridNewGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_LHC_${MODEL}.json \
    --cycles 0 \
    --output \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'interactive' \
    --task-name ${taskname}_hybrid_lhc_finished \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    --X-rtd MINIMIZER_analytic \
    -v1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}_hybrid_lhc.txt
    cd -

elif [[ $MODE == "collect-hybrid-tev" ]]; then

    cd ${defaultdir}/limits_${MODEL}_hybrid_tev/condor

    combineTool.py -M HybridNewGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_LHC_${MODEL}.json \
    --cycles 0 \
    --output \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'interactive' \
    --task-name ${taskname}_hybrid_tev_finished \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    --X-rtd MINIMIZER_analytic \
    -v1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}_hybrid_tev.txt
    cd -

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

elif [[ $MODE == "delete-crashed-jobs" ]]; then
    ############
    # job submission
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    find . -size -4k | grep root > wrong_files.txt
    while read p; do
    echo "Deleting ${p}"
    rm ${p}
    done < wrong_files.txt

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
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_${MODEL}.json \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'condor' \
    --task-name $taskname2 \
    --dry-run \
    --redefineSignalPOI x \
    --setParameterRanges x=0,1 \
    --setParameters r=1 \
    --freezeParameters r -v1 \
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
        title="Classic categorisation 138 fb^{-1} (13 TeV)"
    elif [[ $ANALYSISTYPE == "sm-ml-only" ]]; then
        title="SM categories only 138 fb^{-1} (13 TeV)"
    elif [[ $ANALYSISTYPE == "classic_lowmass" ]]; then
        title="Low-mass categorisation 138 fb^{-1} (13 TeV)"
    elif [[ $ANALYSISTYPE == "with-sm-ml" ]]; then
        title="138 fb^{-1} (13 TeV)"
    else
        title="138 fb^{-1} (13 TeV)"
    fi
    modelname=${MODEL}_13.root
    [[ $OLDFILES == 1 ]] && modelname="${MODEL}_13_old.root"
    for label in "Preliminary" ""; do
        plotLimitGrid.py asymptotic_grid.root \
        --scenario-label="${scenario_label}" \
        --output ${TAG}_${MODEL}_${label} \
        --title-right="${title}" \
        --cms-sub=${label} \
        --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
        --y-range ${y_min},${y_max} \
        --mass_histogram ${sm_like_mass} \
        --mass_histogram_title ${mass_histogram_title} \
        --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${modelname} \
        --x-title "${x_title}" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt
    done
fi
