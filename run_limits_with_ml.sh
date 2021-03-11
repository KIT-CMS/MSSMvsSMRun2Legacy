# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
MODEL=$3
ANALYSISTYPE=$4

if [[ $ANALYSISTYPE == "classic" ]]; then
    analysis="mssm_vs_sm_classic"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_classic"
    fi
else
    analysis="mssm_vs_sm_h125"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_h125"
    fi
fi
# Szenarios from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGMSSMNeutral?redirectedfrom=LHCPhysics.LHCHXSWGMSSMNeutral#Baseline_scenarios
### MSSM scenarios #####
if [[ $MODEL == "mh125" ]]; then
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "mh125_lc" ]]; then
    wsoutput="ws_mh125_lc.root"
    modelfile="13,Run2017,mh125_lc_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "mh125_ls" ]]; then
    wsoutput="ws_mh125_ls.root"
    modelfile="13,Run2017,mh125_ls_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h}^{125}(#tilde{#tau}) scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "mh125_align" ]]; then
    wsoutput="ws_mh125_align.root"
    modelfile="13,Run2017,mh125_align_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h}^{125} alignment scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "HH125" ]]; then
    wsoutput="ws_mHH125.root"
    modelfile="13,Run2017,mHH125_13.root"
    min_mass=150
    max_mass=200
    scenario_label="M_{H}^{125} alignment scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "hm1125_cpv" ]]; then
    #TODO this still need more modifications in the code
    wsoutput="ws_hm1125_cpv.root"
    modelfile="13,Run2017,hm1125_CPV_13.root"
    analysis="mssm_vs_sm_CPV"
    min_mass=130
    max_mass=1500
    scenario_label="M_{h_1}^{125} (CPV) scenario (h,H,A#rightarrow#tau#tau)"
### Negative mu scenarios #####
elif [[ $MODEL == "mh125_muneg_1" ]]; then
    wsoutput="mh125_muneg_1.root"
    modelfile="13,Run2017,mh125_muneg_1_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h}^{125} (#mu = -1 TeV) scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "mh125_muneg_2" ]]; then
    wsoutput="mh125_muneg_2.root"
    modelfile="13,Run2017,mh125_muneg_2_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h}^{125} (#mu = -2 TeV) scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "mh125_muneg_3" ]]; then
    wsoutput="mh125_muneg_3.root"
    modelfile="13,Run2017,mh125_muneg_3_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h}^{125} (#mu = -3 TeV) scenario (h,H,A#rightarrow#tau#tau)"
### EFT scenarios #####
elif [[ $MODEL == "mh125EFT" ]]; then
    wsoutput="mh125EFT.root"
    modelfile="13,Run2017,mh125EFT_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h,#text{EFT}}^{125} scenario (h,H,A#rightarrow#tau#tau)"
elif [[ $MODEL == "mh125EFT_lc" ]]; then
    wsoutput="mh125EFT_lc.root"
    modelfile="13,Run2017,mh125EFT_lc_13.root"
    min_mass=110
    max_mass=3200
    scenario_label="M_{h,#text{EFT}}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
else
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    min_mass=110
    max_mass=3200
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
            --eras 2016,2017,2018 \
            --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt\
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --eras 2016,2017,2018 \
            --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt\
            --variable mt_tot_puppi \
            --parallel 1 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_control_log.txt
    else
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --eras 2016,2017,2018 \
            --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_neuralnet_categories.txt \
            --variable nnscore \
            --sm \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --eras 2016,2017,2018 \
            --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_new_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --eras 2016,2017,2018 \
            --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt \
            --variable mt_tot_puppi \
            --parallel 1 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
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
    --PO modelFile=${modelfile} \
    --PO minTemplateMass=${min_mass} \
    --PO maxTemplateMass=${max_mass} \
    --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v0.root \
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
    --redefineSignalPOI x \
    --setParameters r=1 \
    --freezeParameters r -v1 \
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

elif [[ $MODE == "submit-gc" ]]; then
    ############
    # job submission
    ############
    python scripts/build_gc_job.py \
        --combine-script ${defaultdir}/limits_${MODEL}/condor/condor_${taskname}.sh \
        --workspace ${datacarddir}/combined/cmb/${wsoutput} \
        --workdir /work/sbrommer/workdirs/combine/${taskname} \
        --tag ${taskname} \
        --se-path /storage/gridka-nrg/sbrommer/gc_storage/combine/${TAG}/${taskname}

    ${CMSSW_BASE}/src/grid-control/go.py /work/sbrommer/workdirs/combine/${taskname}/${taskname}.conf -Gc -m 3

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
    --setParameters r=1 \
    --freezeParameters r -v1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/collect_jobs_${MODEL}.txt

    condor_submit condor_${taskname2}.sub
    cp asymptotic_grid.root ..
    cd ${defaultdir}/limits_${MODEL}/

    ############
    # limit plot
    ############
    if [[ $ANALYSISTYPE == "classic" ]]; then
        title="Classic categorisation 137 fb^{-1} (13 TeV)"
    else
        title="ML categorisation 137 fb^{-1} (13 TeV)"
    fi
    plotLimitGrid.py asymptotic_grid.root \
    --scenario-label="${scenario_label}" \
    --output ${TAG}_${MODEL} \
    --title-right="${title}" \
    --cms-sub="Own Work" \
    --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
    --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${MODEL}_13.root \
    --y-range 2.0,60.0 \
    --x-title "m_{A} [GeV]" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt
fi
