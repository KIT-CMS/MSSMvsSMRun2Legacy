# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
if [[ $TAG == "auto" ]]; then
    TAG="cmb_ind"
fi


defaultdir=vlq_analysis/$TAG
analysis="vector_leptoquark"
hSM_treatment="hSM-in-bg"
categorization="classic"
sm_like_hists="bsm"
sub_analysis="betaRd33_0_offdiag0"
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/limits/condor ]] && mkdir -p ${defaultdir}/limits/condor
[[ ! -d ${defaultdir}/limits_ind/condor ]] && mkdir -p ${defaultdir}/limits_ind/condor
defaultdir=$(readlink -f analysis/$TAG)
datacarddir=${defaultdir}/datacards_${analysis}
freeze="MH=200,r_ggH=0.0,r_bbH=0.0"

if [[ $MODE == "initial" ]]; then
    ############
    # morphing
    ############
    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --sub-analysis ${sub_analysis} \
        --hSM-treatment ${hSM_treatment} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_fake_categories.txt \
        --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1" \
        --variable mt_tot_puppi \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_vlq_log.txt

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

    combineTool.py -M T2W -o "ws.root" \
    -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
    -i ${datacarddir}/combined/cmb/ \
    -m 1000.0 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_ind/condor
    combineTool.py -m "500,1000,2000,3000,4000,5000" \
    -M AsymptoticLimits \
    --setParameters gU=0 \
    --redefineSignalPOIs gU \
    -d ${datacarddir}/combined/cmb/ws.root \
    --there -n ".limit" \
    --job-mode condor \
    --dry-run \
    --task-name vlq_full_cmb \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh.txt

elif [[ $MODE == "ws_sge" ]]; then
    ############
    # workspace creation
    ############

    combineTool.py -M T2W -o "ws.root" \
    -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
    -i ${datacarddir}/combined/cmb/ \
    -m 1000 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

    ############
    # job setup creation
    ############
    combineTool.py -m "500,1000,2000,3000,4000,5000" \
    -M AsymptoticLimits \
    --setParameters gU=0 \
    --redefineSignalPOIs gU \
    -d ${datacarddir}/combined/cmb/ws.root \
    --there -n ".limit" \
    --job-mode "SGE" \
    --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
    --task-name vlq_${sub_analysis}_full_cmb \
    --X-rtd MINIMIZER_analytic \
    -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh.txt


elif [[ $MODE == "collect" ]]; then
    combineTool.py -M CollectLimits ${datacarddir}/combined/cmb/higgsCombine.limit*.root \
    --use-dirs \
    -o ${datacarddir}/combined/cmb/vlq_${sub_analysis}.json
 
    plotMSSMLimits.py --cms-sub "Preliminary" \
    --title-right "138 fb^{-1} (13 TeV)" \
    --x-title "M_{U} [TeV]"\
    --y-axis-min 0.0 \
    --y-axis-max 6.0 \
    --show exp ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
    --process "vector_leptoquark" \
    --output vlq_${sub_analysis}_cmb
fi
