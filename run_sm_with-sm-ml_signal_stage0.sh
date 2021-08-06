# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
if [[ $TAG == "auto" ]]; then
    TAG="cmb_sm_analysis_allcats"
fi


defaultdir=analysis/$TAG
analysis="sm"
sub_analysis="no-hSM-in-bg"
categorization="with-sm-ml"
sm_like_hists="sm125"
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/signal_strength/condor ]] && mkdir -p ${defaultdir}/signal_strength/condor
defaultdir=$(readlink -f analysis/$TAG)
datacarddir=${defaultdir}/datacards_${analysis}

if [[ $MODE == "initial" ]]; then
    ############
    # morphing
    ############
    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --sub-analysis ${sub_analysis} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} --sm \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_neuralnet_categories.txt \
        --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1" \
        --variable nnscore \
        --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_sm_cats_log.txt

    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --sub-analysis ${sub_analysis} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_new_categories.txt \
        --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1" \
        --variable mt_tot_puppi \
        --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_mssm_cats_log.txt

    ############
    # combining outputs
    ############
    mkdir -p ${datacarddir}/combined/cmb/
    rsync -av --progress ${datacarddir}/201?/htt_*/* ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
    for ERA in 2016 2017 2018; do
        for CH in et mt tt em; do
            mkdir -p ${datacarddir}/${ERA}/${CH}/
            rsync -av --progress ${datacarddir}/${ERA}/htt_${CH}*/* ${datacarddir}/${ERA}/${CH}/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
        done
        mkdir -p ${datacarddir}/${ERA}/cmb/
        rsync -av --progress ${datacarddir}/${ERA}/htt_*/* ${datacarddir}/${ERA}/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
    done

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############

    combineTool.py -M T2W -o "ws_sm.root" \
    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
    --PO '"map=^.*/(bb|gg)H125.?$:r_ggH[1,0,2]"' \
    --PO '"map=^.*/(W|Z|qq)H125.?$:r_qqH[1,0,2]"' \
    -i ${datacarddir}/combined/cmb/ \
    -m 125.0 --parallel 4 | tee -a ${defaultdir}/logs/workspace_sm.txt

elif [[ $MODE == "run" ]]; then
    ############
    # job setup creation
    ############
    cd ${defaultdir}/signal_strength/condor
    combineTool.py -m "125" \
    -M MultiDimFit \
    --algo singles \
    --robustFit 1 \
    -d ${datacarddir}/combined/cmb/ws_sm.root \
    --there -n ".stage0" \
    --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    --floatOtherPOIs 1 \
    --setParameters r_ggH=1,r_qqH=1 --setParameterRanges CMS_htt_ttbarShape=-10.0,10.0 \
    -v 1 | tee -a ${defaultdir}/logs/run_stage0_sm.txt
fi
