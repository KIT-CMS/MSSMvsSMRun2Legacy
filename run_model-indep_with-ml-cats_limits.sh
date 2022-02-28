# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
if [[ $TAG == "auto" ]]; then
    TAG="cmb_ind_with-sm-ml"
fi


defaultdir=analysis/$TAG
analysis="bsm-model-indep"
hSM_treatment="hSM-in-bg"
categorization="with-sm-ml"
sm_like_hists="bsm"
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
    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --sub-analysis "sm-like-light" \
        --hSM-treatment ${hSM_treatment} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_neuralnet_categories.txt \
        --variable nnscore \
        --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --sm \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --hSM-treatment ${hSM_treatment} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_new_categories.txt \
        --variable mt_tot_puppi \
        --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt

    ############
    # combining outputs
    ############
    mkdir -p ${datacarddir}/combined/cmb/
    rsync -av --progress ${datacarddir}/201?/htt_*/* ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############

    combineTool.py -M T2W -o "ws.root" \
    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
    --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' \
    --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' \
    -i ${datacarddir}/combined/cmb/ \
    -m 125.0 --parallel 4 | tee -a ${defaultdir}/logs/workspace_independent.txt

    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_ind/condor
    combineTool.py -m "60,80,100,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500" \
    -M AsymptoticLimits \
    --rAbsAcc 0 \
    --rRelAcc 0.0005 \
    --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json \
    --setParameters r_ggH=0,r_bbH=0 \
    --redefineSignalPOIs r_bbH \
    -d ${datacarddir}/combined/cmb/ws.root \
    --there -n ".bbH" \
    --job-mode condor \
    --dry-run \
    --task-name bbH_full_cmb \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh.txt

    combineTool.py -m "60,80,100,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500" \
    -M AsymptoticLimits \
    --rAbsAcc 0 \
    --rRelAcc 0.0005 \
    --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json \
    --setParameters r_ggH=0,r_bbH=0 \
    --redefineSignalPOIs r_ggH \
    -d ${datacarddir}/combined/cmb/ws.root \
    --there -n ".ggH" \
    --job-mode condor \
    --dry-run \
    --task-name ggH_full_cmb \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_ggh.txt

elif [[ $MODE == "submit" ]]; then
    ############
    # job submission
    ############
    cd ${defaultdir}/limits_ind/condor
    condor_submit condor_bbH_full_cmb.sub
    condor_submit condor_ggH_full_cmb.sub

    elif [[ $MODE == "submit-local" ]]; then
    ############
    # job submission
    ############
    cp scripts/run_limits_locally.py ${defaultdir}/limits_ind/condor
    cd ${defaultdir}/limits_ind/condor
    python run_limits_locally.py --cores 10 --njobs 31 --taskname condor_bbH_full_cmb.sh
    python run_limits_locally.py --cores 10 --njobs 31 --taskname condor_ggH_full_cmb.sh

elif [[ $MODE == "collect" ]]; then
    for p in gg bb
    do
        combineTool.py -M CollectLimits ${datacarddir}/combined/cmb/higgsCombine.${p}H*.root \
        --use-dirs \
        -o ${datacarddir}/combined/cmb/mssm_${p}H_cmb.json

        plotMSSMLimits.py --cms-sub "Preliminary" \
        --title-right "138 fb^{-1} (13 TeV)" \
        --process "${p}#phi" \
        --y-axis-min 0.0001 \
        --y-axis-max 1000.0 \
        --show exp,obs ${datacarddir}/combined/cmb/mssm_${p}H_$cmb.json  \
        --output mssm_model-independent_${p}H_cmb \
        --logx \
        --logy
    done
fi
