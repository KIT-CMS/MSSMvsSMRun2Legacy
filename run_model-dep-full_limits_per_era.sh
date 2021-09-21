# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
ERA=$3
ANALYSISTYPE=$4
HSMTREATMENT=$5

if [[ $ANALYSISTYPE == "classic" ]]; then
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="classic"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_${ERA}_classic"
    fi
else
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="with-sm-ml"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_${ERA}_with_ml"
    fi
fi

if [[ $ERA == "2016" ]]; then
    LUMI="36.3 fb^{-1} (2016, 13 TeV)"
elif [[ $ERA == "2017" ]]; then
    LUMI="41.5 fb^{-1} (2017, 13 TeV)"
elif [[ $ERA == "2018" ]]; then
    LUMI="59.7 fb^{-1} (2018, 13 TeV)"
fi

defaultdir=analysis/$TAG
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/limits/condor ]] && mkdir -p ${defaultdir}/limits/condor
[[ ! -d ${defaultdir}/limits_ind/condor ]] && mkdir -p ${defaultdir}/limits_ind/condor
defaultdir=$(readlink -f analysis/$TAG)
datacarddir=${defaultdir}/datacards_${analysis}
taskname="${analysis}_${TAG}_1"
taskname2="${analysis}_${TAG}_2"

if [[ $MODE == "initial" ]]; then
    ############
    # morphing
    ############
    if [[ $ANALYSISTYPE == "classic" ]]; then
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment ${HSMTREATMENT} \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
    else
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment ${HSMTREATMENT} \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_neuralnet_categories.txt \
            --variable nnscore \
            --sm \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment ${HSMTREATMENT} \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
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
    mkdir -p ${datacarddir}/${ERA}/cmb/
    rsync -av --progress ${datacarddir}/${ERA}/htt_*/* ${datacarddir}/${ERA}/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards_${ERA}.txt

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############
    combineTool.py -M T2W -o ws_mh125.root \
    -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM \
    --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
    --PO replace-with-SM125=${replace_with_sm125} \
    --PO modelFile=13,Run2017,mh125_13.root \
    --PO minTemplateMass=60 \
    --PO maxTemplateMass=3500 \
    --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
    -i ${datacarddir}/$ERA/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace.txt
    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits/condor
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json \
    -d ${datacarddir}/$ERA/cmb/ws_mh125.root \
    --job-mode 'condor' \
    --task-name $taskname \
    --dry-run \
    --redefineSignalPOI x \
    --setParameters r=1 \
    --freezeParameters r -v1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 &> ${defaultdir}/logs/job_setup.txt

elif [[ $MODE == "submit" ]]; then
    ############
    # job submission
    ############
    cd ${defaultdir}/limits/condor
    condor_submit condor_${taskname}.sub

elif [[ $MODE == "submit-local" ]]; then
    ############
    # job submission
    ############
    cp scripts/run_limits_locally.py ${defaultdir}/limits/condor
    cd ${defaultdir}/limits/condor
    python run_limits_locally.py --cores 12 --taskname condor_${taskname}.sh

elif [[ $MODE == "collect" ]]; then
    ############
    # job collection
    ############
    cd ${defaultdir}/limits/condor/
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json \
    -d ${datacarddir}/$ERA/cmb/ws_mh125.root \
    --job-mode 'condor' \
    --task-name $taskname2 \
    --dry-run \
    --redefineSignalPOI x \
    --setParameters r=1 \
    --freezeParameters r -v1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/collect_jobs.txt

    condor_submit condor_${taskname2}.sub
    cp asymptotic_grid.root ..
    cd ${defaultdir}/limits/

    ############
    # limit plot
    ############
    plotLimitGrid.py asymptotic_grid.root \
    --scenario-label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)" \
    --output ${TAG}_${ERA}_cmb \
    --title-right="${LUMI}" \
    --cms-sub="Preliminary" \
    --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
    --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/mh125_13.root \
    --y-range 2.0,60.0 \
    --x-title "m_{A} [GeV]" 2>&1 | tee -a ${defaultdir}/logs/plot_grid.txt
fi
