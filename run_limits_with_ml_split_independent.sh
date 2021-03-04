# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
CHANNEL=$3
ERA=$4
if [[ $TAG == "auto" ]]; then
    TAG="${CHANNEL}_${ERA}_h125"
fi

if [[ $ERA == "2016" ]]; then
    LUMI="35.9 fb^{-1} (2016, 13 TeV)"
elif [[ $ERA == "2017" ]]; then
    LUMI="41.5 fb^{-1} (2017, 13 TeV)"
elif [[ $ERA == "2018" ]]; then
    LUMI="59.7 fb^{-1} (2018, 13 TeV)"
fi

defaultdir=analysis/$TAG
analysis="mssm"
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
        --eras $ERA \
        --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/by_channel/sm_neuralnet_categories_$CHANNEL.txt \
        --variable nnscore \
        --sm_gg_fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v0.root \
        --sm \
        --parallel 1 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

    morph_parallel.py --output ${defaultdir}/datacards \
        --analysis ${analysis} \
        --eras $ERA \
        --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/by_channel/mssm_signal_categories_$CHANNEL.txt \
        --variable mt_tot_puppi \
        --sm_gg_fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v0.root \
        --parallel 6 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt

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

    combineTool.py -M T2W -o "ws.root" \
    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
    --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' \
    --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' \
    -i ${datacarddir}/${ERA}/${CHANNEL} \
    -m 110 --parallel 4 | tee -a ${defaultdir}/logs/workspace_independent_${ERA}_${CHANNEL}.txt


    # # Extract prefit shapes.
    # prefit_postfit_shapes_parallel.py --datacard_pattern "output_mssm/201?/htt_*/combined.txt.cmb" \
    #                                 --workspace_name ws.root \
    #                                 --freeze_arguments "--freeze $FREEZE" \
    #                                 --output_name prefit_shapes_${FREEZE}.root
    #                                 --parallel $CORES

    # hadd -f prefit_shapes_${FREEZE}.root output_mssm/201?/htt_*/prefit_shapes_${FREEZE}.root

    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_ind/condor
    combineTool.py -m "110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" \
    -M AsymptoticLimits \
    --rAbsAcc 0 \
    --rRelAcc 0.0005 \
    --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json \
    --setParameters r_ggH=0,r_bbH=0 \
    --redefineSignalPOIs r_bbH \
    -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
    --there -n ".bbH" \
    --job-mode condor \
    --dry-run \
    --task-name bbH_full_${ERA}_${CHANNEL} \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh_${ERA}_${CHANNEL}.txt


    combineTool.py -m "110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" \
    -M AsymptoticLimits \
    --rAbsAcc 0 \
    --rRelAcc 0.0005 \
    --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json \
    --setParameters r_ggH=0,r_bbH=0 \
    --redefineSignalPOIs r_ggH \
    -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
    --there -n ".ggH" \
    --job-mode condor \
    --dry-run \
    --task-name ggH_full_${ERA}_${CHANNEL} \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerStrategy 0 \
    --cminDefaultMinimizerTolerance 0.01 \
    -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_ggh_${ERA}_${CHANNEL}.txt

elif [[ $MODE == "submit" ]]; then
    ############
    # job submission
    ############
    cd ${defaultdir}/limits_ind/condor
    condor_submit condor_bbH_full_${ERA}_${CHANNEL}.sub
    condor_submit condor_ggH_full_${ERA}_${CHANNEL}.sub

    elif [[ $MODE == "submit-local" ]]; then
    ############
    # job submission
    ############
    cp scripts/run_limits_locally.py ${defaultdir}/limits_ind/condor
    cd ${defaultdir}/limits_ind/condor
    python run_limits_locally.py --cores 27 --njobs 27 --taskname condor_bbH_full_${ERA}_${CHANNEL}.sh | tee -a ${defaultdir}/logs/fit_bbh_${ERA}_${CHANNEL}.txt
    python run_limits_locally.py --cores 27 --njobs 27 --taskname condor_ggH_full_${ERA}_${CHANNEL}.sh | tee -a ${defaultdir}/logs/fit_ggh_${ERA}_${CHANNEL}.txt

elif [[ $MODE == "collect" ]]; then
    for p in gg bb
    do
        combineTool.py -M CollectLimits ${datacarddir}/${ERA}/${CHANNEL}/higgsCombine.${p}H*.root \
            --use-dirs \
            -o ${datacarddir}/${ERA}/${CHANNEL}/mssm_${p}H_${ERA}_${CHANNEL}.json
    done

elif [[ $MODE == "plot" ]]; then
    for p in gg bb
    do
        plotMSSMLimits.py --cms-sub "Own Work" \
            --title-right="${CHANNEL} - ${LUMI}" \
            --process "${p}#phi" \
            --y-axis-min 0.0001 \
            --y-axis-max 1000.0 \
            --show exp,obs ${datacarddir}/${ERA}/${CHANNEL}/mssm_${p}H_${ERA}_${CHANNEL}_${CHANNEL}.json  \
            --output mssm_model-independent_${p}H_${ERA}_${CHANNEL} \
            --logx \
            --logy
        mv mssm_model-independent_${p}H_${ERA}_${CHANNEL}.{png,pdf} ${defaultdir}/limits_ind
    done

elif [[ $MODE == "prefit" ]]; then
    FREEZE="MH=450,r_ggH=0.002,r_bbH=0.002"
    prefit_postfit_shapes_parallel.py --datacard_pattern "analysis/mt_16_debug_independent/datacards_mssm/$ERA/$CHANNEL/combined.txt.cmb" \
                                  --workspace_name ws.root \
                                  --freeze_arguments "--freeze $FREEZE" \
                                  --output_name prefit_shapes_${FREEZE}.root \
                                  --parallel 1

    hadd -f prefit_shapes_${FREEZE}.root output_mssm/201?/htt_*/prefit_shapes_${FREEZE}.root

elif [[ $MODE == "impacts" ]]; then
    combineTool.py -M Impacts -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
        -m 125 \
        --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
        --doInitialFit --robustFit 1 \
        -t -1 --setParameters r_ggH=0,r_bbH=0 \
        --parallel 16

    combineTool.py -M Impacts -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
        -m 125 \
        --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
        --robustFit 1 --doFits \
        -t -1 --setParameters r_ggH=0,r_bbH=0 \
        --parallel 16

    combineTool.py -M Impacts \
        -d ${datacarddir}/${ERA}/${CHANNEL}/ws.root \
        -m 125 \
        -o ${ERA}_${CHANNEL}_impacts.json

    plotImpacts.py -i ${ERA}_${CHANNEL}_impacts.json -o ${ERA}_${CHANNEL}_impacts
fi
