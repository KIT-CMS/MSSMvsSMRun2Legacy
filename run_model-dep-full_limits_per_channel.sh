# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
CHANNEL=$3
ERA=$4
MODEL=$5
HSMTREATMENT=$6

analysis="bsm-model-dep-full"
sm_like_hists="sm125"
replace_with_sm125=1
categorization="with-sm-ml"
if [[ $TAG == "auto" ]]; then
    TAG="${ERA}_${CHANNEL}_with_ml"
fi
echo $TAG
if [[ $ERA == "2016" ]]; then
    LUMI="36.3 fb^{-1} (2016, 13 TeV)"
elif [[ $ERA == "2017" ]]; then
    LUMI="41.5 fb^{-1} (2017, 13 TeV)"
elif [[ $ERA == "2018" ]]; then
    LUMI="59.7 fb^{-1} (2018, 13 TeV)"
fi
# Szenarios from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGMSSMNeutral?redirectedfrom=LHCPhysics.LHCHXSWGMSSMNeutral#Baseline_scenarios
### MSSM scenarios #####
if [[ $MODEL == "mh125" ]]; then
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    scenario_label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
elif [[ $MODEL == "mh125_lc" ]]; then
    wsoutput="ws_mh125_lc.root"
    modelfile="13,Run2017,mh125_lc_13.root"
    scenario_label="M_{h}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
elif [[ $MODEL == "mh125_ls" ]]; then
    wsoutput="ws_mh125_ls.root"
    modelfile="13,Run2017,mh125_ls_13.root"
    scenario_label="M_{h}^{125}(#tilde{#tau}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
elif [[ $MODEL == "mh125_align" ]]; then
    wsoutput="ws_mh125_align.root"
    modelfile="13,Run2017,mh125_align_13.root"
    scenario_label="M_{h}^{125} alignment scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
elif [[ $MODEL == "mHH125" ]]; then
    wsoutput="ws_mHH125.root"
    modelfile="13,Run2017,mHH125_13.root"
    scenario_label="M_{H}^{125} alignment scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-heavy"
elif [[ $MODEL == "mh1125_CPV" ]]; then
    wsoutput="ws_mh1125_cpv.root"
    modelfile="13,Run2017,mh1125_CPV_13.root"
    scenario_label="M_{h_1}^{125} (CPV) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="cpv"
### Negative mu scenarios #####
elif [[ $MODEL == "mh125_muneg_1" ]]; then
    wsoutput="mh125_muneg_1.root"
    modelfile="13,Run2017,mh125_muneg_1_13.root"
    scenario_label="M_{h}^{125} (#mu = -1 TeV) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
elif [[ $MODEL == "mh125_muneg_2" ]]; then
    wsoutput="mh125_muneg_2.root"
    modelfile="13,Run2017,mh125_muneg_2_13.root"
    scenario_label="M_{h}^{125} (#mu = -2 TeV) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
elif [[ $MODEL == "mh125_muneg_3" ]]; then
    wsoutput="mh125_muneg_3.root"
    modelfile="13,Run2017,mh125_muneg_3_13.root"
    scenario_label="M_{h}^{125} (#mu = -3 TeV) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
### EFT scenarios #####
elif [[ $MODEL == "mh125EFT" ]]; then
    wsoutput="mh125EFT.root"
    modelfile="13,Run2017,mh125EFT_13.root"
    scenario_label="M_{h,#text{EFT}}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
elif [[ $MODEL == "mh125EFT_lc" ]]; then
    wsoutput="mh125EFT_lc.root"
    modelfile="13,Run2017,mh125EFT_lc_13.root"
    scenario_label="M_{h,#text{EFT}}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
else
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    sub_analysis="sm-like-light"
fi
defaultdir=analysis/$TAG
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
defaultdir=$(readlink -f analysis/$TAG)
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/prefitplots ]] && mkdir -p ${defaultdir}/prefitplots
[[ ! -d ${defaultdir}/limits_${MODEL}/condor ]] && mkdir -p ${defaultdir}/limits_${MODEL}/condor

datacarddir=${defaultdir}/datacards_${analysis}
taskname="${analysis}_${TAG}_1"
taskname2="${analysis}_${TAG}_2"

if [[ $MODE == "initial" ]]; then
    ############
    # morphing
    ############
    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --sub-analysis ${sub_analysis} \
        --hSM-treatment ${HSMTREATMENT} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} \
        --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/by_channel/sm_neuralnet_categories_$CHANNEL.txt \
        --variable nnscore \
        --sm \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --sub-analysis ${sub_analysis} \
        --hSM-treatment ${HSMTREATMENT} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} \
        --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/by_channel/mssm_signal_categories_$CHANNEL.txt \
        --variable mt_tot_puppi \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt

    ############
    # combining outputs
    ############
    mkdir -p ${datacarddir}/combined/cmb/
    rsync -av --progress ${datacarddir}/201?/htt_*/* ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
    mkdir -p ${datacarddir}/${ERA}/cmb/

    rsync -av --progress ${datacarddir}/${ERA}/htt_*/* ${datacarddir}/${ERA}/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards_${ERA}.txt
    mkdir -p ${datacarddir}/${ERA}/${CHANNEL}/
    rsync -av --progress ${datacarddir}/${ERA}/htt_${CHANNEL}*/* ${datacarddir}/${ERA}/${CHANNEL}/
    rsync -av --progress ${datacarddir}/${ERA}/htt_${CHANNEL}* ${datacarddir}/${ERA}/${CHANNEL}/
    rsync -av --progress ${datacarddir}/${ERA}/restore_binning ${datacarddir}/${ERA}/${CHANNEL}/restore_binning

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############
    combineTool.py -M T2W -o ${wsoutput} \
    -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM \
    --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
    --PO replace-with-SM125=${replace_with_sm125} \
    --PO hSM-treatment=$HSMTREATMENT \
    --PO modelFile=${modelfile} \
    --PO minTemplateMass=60 \
    --PO maxTemplateMass=3500 \
    --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
    -i ${datacarddir}/${ERA}/${CHANNEL}/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt

    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_${MODEL}/condor
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_${MODEL}.json \
    -d ${datacarddir}/$ERA/$CHANNEL/${wsoutput}  \
    --job-mode 'condor' \
    --task-name $taskname \
    --dry-run \
    --redefineSignalPOI x \
    --setParameters r=1 \
    --freezeParameters r -v1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 &> ${defaultdir}/logs/job_setup_${MODEL}.txt

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
    python run_limits_locally.py --cores 32 --taskname condor_${taskname}.sh

elif [[ $MODE == "hybrid-gc" ]]; then

    mkdir -p ${defaultdir}/limits_${MODEL}_hybrid/condor
    cd ${defaultdir}/limits_${MODEL}_hybrid/condor

    combineTool.py -M HybridNewGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_${MODEL}.json \
    --cycles 2 \
    -d ${datacarddir}/$ERA/$CHANNEL/${wsoutput}  \
    --job-mode 'condor' \
    --task-name ${taskname}_hybrid \
    --dry-run | tee -a ${defaultdir}/logs/job_setup_${MODEL}_hybrid.txt
    ############
    # job submission
    ############
    cd -

elif [[ $MODE == "check-hybrid-gc" ]]; then

    cd ${defaultdir}/limits_${MODEL}_hybrid/condor

    combineTool.py -M HybridNewGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_hybrid_grid_${MODEL}.json \
    --cycles 0 \
    -d ${datacarddir}/$ERA/$CHANNEL/${wsoutput} | tee -a ${defaultdir}/logs/job_setup_${MODEL}_hybrid.txt
    ############
    # job submission
    ############
    cd -

elif [[ $MODE == "submit-gc" ]]; then
    ############
    # job submission
    ############
    python scripts/build_gc_job.py \
        --combine-script ${defaultdir}/limits_${MODEL}/condor/condor_${taskname}.sh \
        --workspace ${datacarddir}}/$ERA/$CHANNEL/${wsoutput} \
        --workdir /work/sbrommer/workdirs/combine/${taskname} \
        --tag ${taskname} \
        --se-path /storage/gridka-nrg/sbrommer/gc_storage/combine/${TAG}/${taskname}

    ${CMSSW_BASE}/src/grid-control/go.py /work/sbrommer/workdirs/combine/${taskname}/${taskname}.conf -Gc -m 3

elif [[ $MODE == "copy-results-gc" ]]; then
    ############
    # job submission
    ############
    rsync -avhP /storage/gridka-nrg/sbrommer/gc_storage/combine/${TAG}/${taskname}/output/ ${defaultdir}/limits_${MODEL}/condor


elif [[ $MODE == "collect" ]]; then
    ############
    # job collection
    ############
    cd ${defaultdir}/limits_${MODEL}/condor/
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json \
    -d ${datacarddir}/$ERA/$CHANNEL/${wsoutput}  \
    --job-mode 'condor' \
    --task-name $taskname2 \
    --dry-run \
    --redefineSignalPOI x \
    --setParameters r=1 \
    --freezeParameters r -v1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/collect_jobs_${MODEL}.txt

    cp asymptotic_grid.root ..
    cd ${defaultdir}/limits_${MODEL}/

    ############
    # limit plot
    ############
    plotLimitGrid.py asymptotic_grid.root \
    --scenario-label="${scenario_label}" \
    --output ${TAG}_${ERA}_${CHANNEL} \
    --title-right="${CHANNEL} - ${LUMI}" \
    --cms-sub="Preliminary" \
    --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
    --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${MODEL}_13.root \
    --y-range 2.0,60.0 \
    --x-title "m_{A} [GeV]" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt

elif [[ $MODE == "plot-prefit" ]]; then
    combineTool.py -M T2W -o ${wsoutput} \
        -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM \
        --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
        --PO modelFile=${modelfile} \
        --PO minTemplateMass=${min_mass} \
        --PO maxTemplateMass=${max_mass} \
        --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        -i ${datacarddir}/201?/*/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt
    source utils/setup_python.sh
    python scripts/prefit_shapes_v2.py --basedir ${datacarddir} --workspacename ${wsoutput}
    hadd -f  ${defaultdir}/prefitplots/prefitshapes.root ${datacarddir}/prefitshapes/htt_*
        for OPTION in "" "--png"
        do
            ./plotting/plot_shapes_mssm.py -i ${defaultdir}/prefitplots/prefitshapes.root -c ${CHANNEL} \
                                        -e $ERA $OPTION --fake-factor --embedding --normalize-by-bin-width \
                                        -o ${defaultdir}/prefitplots --use-sm --linear --blinded
        done
fi
