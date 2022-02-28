# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
MODEL=$3
ANALYSISTYPE=$4
GRIDUSER=$5

if [[ $ANALYSISTYPE == "classic" ]]; then
    analysis="bsm-model-dep-additional"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="classic"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_classic"
    fi
else
    analysis="bsm-model-dep-additional"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="with-sm-ml"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_with_ml"
    fi
fi
# Szenarios from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGMSSMNeutral?redirectedfrom=LHCPhysics.LHCHXSWGMSSMNeutral#Baseline_scenarios
### MSSM scenarios #####
if [[ $MODEL == "mh125" ]]; then
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_lc" ]]; then
    wsoutput="ws_mh125_lc.root"
    modelfile="13,Run2017,mh125_lc_13.root"
    scenario_label="M_{h}^{125}(#tilde{#chi}) scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_ls" ]]; then
    wsoutput="ws_mh125_ls.root"
    modelfile="13,Run2017,mh125_ls_13.root"
    scenario_label="M_{h}^{125}(#tilde{#tau}) scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_align" ]]; then
    wsoutput="ws_mh125_align.root"
    modelfile="13,Run2017,mh125_align_13.root"
    scenario_label="M_{h}^{125} alignment scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=20.0
elif [[ $MODEL == "mHH125" ]]; then
    wsoutput="ws_mHH125.root"
    modelfile="13,Run2017,mHH125_13.root"
    scenario_label="M_{H}^{125} alignment scenario (h,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-heavy"
    sm_like_mass="m_H"
    x_title='m_{H^{#plus}} [GeV]'
    mass_histogram_title="m_{H}"
    y_min=5.0
    y_max=6.0
elif [[ $MODEL == "mh1125_CPV" ]]; then
    wsoutput="ws_mh1125_cpv.root"
    modelfile="13,Run2017,mh1125_CPV_13.root"
    scenario_label="M_{h_{1}}^{125} (CPV) scenario (^{}h_{2},^{}h_{3}#rightarrow#tau#tau)"
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
    scenario_label="M_{h}^{125} (#mu = -1 TeV) scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=56.0
elif [[ $MODEL == "mh125_muneg_2" ]]; then
    wsoutput="mh125_muneg_2.root"
    modelfile="13,Run2017,mh125_muneg_2_13.root"
    scenario_label="M_{h}^{125} (#mu = -2 TeV) scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=30.0
elif [[ $MODEL == "mh125_muneg_3" ]]; then
    wsoutput="mh125_muneg_3.root"
    modelfile="13,Run2017,mh125_muneg_3_13.root"
    scenario_label="M_{h}^{125} (#mu = -3 TeV) scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=20.0
### EFT scenarios #####
elif [[ $MODEL == "mh125EFT" ]]; then
    wsoutput="mh125EFT.root"
    modelfile="13,Run2017,mh125EFT_13.root"
    scenario_label="M_{h,EFT}^{125} scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=10.0
elif [[ $MODEL == "mh125EFT_lc" ]]; then
    wsoutput="mh125EFT_lc.root"
    modelfile="13,Run2017,mh125EFT_lc_13.root"
    scenario_label="M_{h,EFT}^{125}(#tilde{#chi}) scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=10.0
else
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    scenario_label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} [GeV]'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
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
            --hSM-treatment "hSM-in-bg" \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1" \
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
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1" \
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
    -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM \
    --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
    --PO replace-with-SM125=${replace_with_sm125} \
    --PO modelFile=${modelfile} \
    --PO minTemplateMass=60 \
    --PO maxTemplateMass=3500 \
    --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
    --PO sm-predictions=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_predictions_13TeV.json \
    -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt

    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
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

elif [[ $MODE == "hybrid-gc" ]]; then

    mkdir -p ${defaultdir}/limits_${MODEL}_hybrid/condor
    cd ${defaultdir}/limits_${MODEL}_hybrid/condor

    combineTool.py -M HybridNewGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_${MODEL}.json \
    --cycles 2 \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --job-mode 'condor' \
    --task-name ${taskname}_hybrid \
    --dry-run | tee -a ${defaultdir}/logs/job_setup_${MODEL}_hybrid.txt
    ############
    # job submission
    ############
    cd -
    # python scripts/build_gc_job.py \
    #     --combine-script ${defaultdir}/limits_${MODEL}_hybrid/condor/condor_${taskname}_hybrid.sh \
    #     --workspace ${datacarddir}/combined/cmb/${wsoutput} \
    #     --workdir /work/sbrommer/workdirs/combine/${taskname}_hybrid \
    #     --tag ${taskname}_hybrid \
    #     --se-path /storage/gridka-nrg/sbrommer/gc_storage/combine/${TAG}/${taskname}_hybrid

    # ${CMSSW_BASE}/src/grid-control/go.py /work/sbrommer/workdirs/combine/${taskname}_hybrid/${taskname}_hybrid.conf -Gc -m 3

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
    --redefineSignalPOI r \
    --setParameterRanges r=0,1 \
    --setParameters r=1,x=1 \
    --freezeParameters x -v1 \
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
    else
        title="138 fb^{-1} (13 TeV)"
    fi
    plotLimitGrid.py asymptotic_grid.root \
    --scenario-label="${scenario_label}" \
    --output ${TAG}_${MODEL} \
    --title-right="${title}" \
    --cms-sub="Preliminary" \
    --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
    --y-range ${y_min},${y_max} \
    --mass_histogram ${sm_like_mass} \
    --mass_histogram_title ${mass_histogram_title} \
    --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${MODEL}_13.root \
    --x-title "${x_title}" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt
fi
