# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
ANALYSISTYPE=$3

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

defaultdir=analysis/$TAG
analysis="mssm_vs_sm_h125"
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

elif [[ $MODE == "ws-dependent" ]]; then
    ############
    # workspace creation
    ############
    combineTool.py -M T2W -o ws_mh125.root \
    -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM \
    --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
    --PO modelFile=13,Run2017,mh125_13.root \
    --PO minTemplateMass=110.0 \
    --PO maxTemplateMass=3200.0 \
    --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v0.root \
    -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace.txt
    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits/condor
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json \
    -d ${datacarddir}/combined/cmb/ws_mh125.root \
    --job-mode 'condor' \
    --task-name $taskname \
    --dry-run \
    --redefineSignalPOI x \
    --setParameters r=1 \
    --freezeParameters r -v1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/job_setup.txt

elif [[ $MODE == "submit-dependent" ]]; then
    ############
    # job submission
    ############
    cd ${defaultdir}/limits/condor
    condor_submit condor_${taskname}.sub

elif [[ $MODE == "submit-dependent-local" ]]; then
    ############
    # job submission
    ############
    cp scripts/run_limits_locally.py ${defaultdir}/limits/condor
    cd ${defaultdir}/limits/condor
    python run_limits_locally.py --cores 20 --taskname condor_${taskname}.sh

elif [[ $MODE == "collect-dependent" ]]; then
    ############
    # job collection
    ############
    cd ${defaultdir}/limits/condor/
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json \
    -d ${datacarddir}/combined/cmb/ws_mh125.root \
    --job-mode 'condor' \
    --task-name $taskname2 \
    --dry-run \
    --redefineSignalPOI x \
    --setParameters r=1 \
    --freezeParameters r -v1 \
    --cminDefaultMinimizerStrategy 0 \
    --X-rtd MINIMIZER_analytic \
    --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/collect_jobs.txt

    condor_submit ${taskname2}.sub
    cp asymptotic_grid.root ..
    cd ${defaultdir}/limits/

    ############
    # limit plot
    ############
    plotLimitGrid.py asymptotic_grid.root \
    --scenario-label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)" \
    --output ${TAG} \
    --title-right="137 fb^{-1} (13 TeV)" \
    --cms-sub="Own Work" \
    --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
    --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/mh125_13.root \
    --y-range 2.0,60.0 \
    --x-title "m_{A} [GeV]" 2>&1 | tee -a ${defaultdir}/logs/plot_grid.txt

elif [[ $MODE == "ws-independent" ]]; then
    ############
    # workspace creation
    ############
    combineTool.py -M T2W -o "ws.root" \
    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
    --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' \
    --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' \
    -i ${datacarddir}/{2016,2017,2018,combined}/cmb \
    -m 110 --parallel 4 | tee -a ${defaultdir}/logs/workspace_independent.txt
    for era in 2016 2017 2018;
    do
        for channel in "et" "mt" "tt";
        do
            combineTool.py -M T2W -o "ws.root" \
            -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
            --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' \
            --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' \
            -i ${datacarddir}/${era}/${channel} \
            -m 110 --parallel 4 | tee -a ${defaultdir}/logs/workspace_independent_${era}_${channel}.txt
        done
    done
    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_ind/condor
    for era in "2016" "2017" "2018" "combined";
    do
        for channel in "et" "mt" "tt" "cmb";
        do
            [[ ! -d ${defaultdir}/limits_ind/condor/${era}/${channel} ]] && mkdir -p ${defaultdir}/limits_ind/condor/${era}/${channel}
            cd ${defaultdir}/limits_ind/condor/${era}/${channel}
            combineTool.py -m "110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" \
            -M AsymptoticLimits \
            --rAbsAcc 0 \
            --rRelAcc 0.0005 \
            --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json \
            --setParameters r_ggH=0,r_bbH=0 \
            --redefineSignalPOIs r_bbH \
            -d ${datacarddir}/${era}/${channel}/ws.root \
            --there -n ".bbH" \
            --job-mode condor \
            --dry-run \
            --task-name bbH_full_${era}_${channel} \
            --X-rtd MINIMIZER_analytic \
            --cminDefaultMinimizerStrategy 0 \
            --cminDefaultMinimizerTolerance 0.01 \
            -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh_${era}_${channel}.txt

            combineTool.py -m "110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" \
            -M AsymptoticLimits \
            --rAbsAcc 0 \
            --rRelAcc 0.0005 \
            --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json \
            --setParameters r_ggH=0,r_bbH=0 \
            --redefineSignalPOIs r_ggH \
            -d ${datacarddir}/${era}/${channel}/ws.root \
            --there -n ".ggH" \
            --job-mode condor \
            --dry-run \
            --task-name ggH_full_${era}_${channel} \
            --X-rtd MINIMIZER_analytic \
            --cminDefaultMinimizerStrategy 0 \
            --cminDefaultMinimizerTolerance 0.01 \
            -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_ggh_${era}_${channel}.txt
        done
    done
elif [[ $MODE == "submit-independent" ]]; then
    ############
    # job submission
    ############

    for era in "2018";
    do
        for channel in "et" "mt" "tt";
        do
            cd ${defaultdir}/limits_ind/condor/${era}/${channel}
            condor_submit condor_bbH_full_${era}_${channel}.sub
            condor_submit condor_ggH_full_${era}_${channel}.sub
        done
    done

elif [[ $MODE == "collect-independent" ]]; then
    for era in "2018";
    do
        for channel in "et" "mt" "tt";
        do
            for p in gg bb
            do
                combineTool.py -M CollectLimits ${datacarddir}/${era}/${channel}/higgsCombine.${p}H*.root \
                --use-dirs \
                -o ${datacarddir}/${era}/${channel}/mssm_${p}H_${era}_${channel}.json

                plotMSSMLimits.py --cms-sub "Own Work" \
                --title-right "137 fb^{-1} (13 TeV)" \
                --process "${p}#phi" \
                --y-axis-min 0.0001 \
                --y-axis-max 1000.0 \
                --show exp,obs ${datacarddir}/${era}/${channel}/mssm_${p}H_${era}_${channel}_${channel}.json  \
                --output ${datacarddir}/limits_ind/mssm_model-independent_${p}H_${era}_${channel} \
                --logx \
                --logy
            done
        done
    done
fi
