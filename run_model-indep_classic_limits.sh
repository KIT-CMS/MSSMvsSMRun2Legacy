#!/bin/bash
# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
if [[ $TAG == "auto" ]]; then
    TAG="cmb_ind"
fi


defaultdir=analysis/$TAG
analysis="bsm-model-indep"
hSM_treatment="hSM-in-bg"
categorization="classic"
sm_like_hists="bsm"
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
        --sub-analysis "sm-like-light" \
        --hSM-treatment ${hSM_treatment} \
        --categorization ${categorization} \
        --sm-like-hists ${sm_like_hists} \
        --eras 2016,2017,2018 \
        --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt \
        --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1" \
        --variable mt_tot_puppi \
        --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt

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

elif [[ $MODE == "ws-gof" ]]; then
    for CH in et mt tt em; do
        rsync -av --progress ${datacarddir}/201?/htt_${CH}_*/* ${datacarddir}/combined/${CH}/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
    done
    # Copy ttbar control region to every channel to include it in the workspaces
    for ERA in 2016 2017 2018; do
        for CH in et mt tt; do
            rsync -av --progress ${datacarddir}/${ERA}/htt_em_2_${ERA}/* ${datacarddir}/${ERA}/${CH}/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
            rsync -av --progress ${datacarddir}/${ERA}/htt_em_2_${ERA}/* ${datacarddir}/combined/${CH}/ 2>&1 | tee -a ${defaultdir}/logs/copy_datacards.txt
        done
    done
    ############
    # workspace creation for GoF tests
    ############

    combineTool.py -M T2W -o "ws-gof.root" \
    -i ${datacarddir}/*/{et,mt,tt,em,cmb}/ \
    --channel-masks \
    -m 125.0 --parallel 8 | tee -a ${defaultdir}/logs/workspace_gof_independent.txt

elif [[ $MODE == "ws-plot" ]]; then
    ###############
    # workspace production for plots
    ###############
    combineTool.py -M T2W -o "ws.root" \
    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
    --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[-50,200]"' \
    --PO '"map=^.*/bbh$:r_bbH[-50,0,200]"' \
    --X-allow-no-signal \
    -i ${datacarddir}/201?/htt_*/ \
    -m 125.0 \
    --parallel 8 | tee -a ${defaultdir}/logs/workspace_plots_independent.txt

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
        --show exp,obs ${datacarddir}/combined/cmb/mssm_${p}H_cmb_cmb.json \
        --output mssm_model-independent_${p}H_cmb \
        --logx \
        --logy
    done

elif [[ $MODE == "prefit-plots" ]]; then
    #####################
    # Extract prefit shapes.
    #####################
    for era in 2016 2017 2018; do
        prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/${era}/htt_em_2_*/combined.txt.cmb" \
                                          --workspace_name ws.root \
                                          --output_name prefit_shapes_${freeze}.root \
                                          --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-combined-${freeze}.log
        prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/${era}/htt_*_3*/combined.txt.cmb" \
                                          --workspace_name ws.root \
                                          --freeze_arguments "--freeze ${freeze}" \
                                          --output_name prefit_shapes_${freeze}.root \
                                          --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-combined-${freeze}.log
    done
    hadd -f ${datacarddir}/combined/cmb/prefit_shapes_${freeze}.root ${datacarddir}/201?/htt_*/prefit_shapes_${freeze}.root | tee -a ${defaultdir}/logs/extract_model_independent_shapes-combined-${freeze}.log

    for era in 2016 2017 2018; do
        bash plotting/plot_shapes_mssm_model_independent.sh \
            ${era} \
            "${datacarddir}/combined/cmb/prefit_shapes_${freeze}.root" \
            "${datacarddir}/plots/prefit_shapes_$(echo ${freeze} | sed 's/=//g; s/\./p/g')/" \
            et,mt,tt,em \
            $(echo $freeze | cut -d, -f1 | cut -d= -f2) \
            $(echo $freeze | cut -d, -f2 | cut -d= -f2)
    done

elif [[ $MODE == "fit-for-plots" ]]; then
    combineTool.py -M FitDiagnostics \
        -d ${datacarddir}/combined/cmb/ws.root \
        -m 200 \
        --setParameters r_ggH=0,r_bbH=0 --setParameterRange r_ggH=-0.00001,0.00001:r_bbH=-2,5 \
        --redefineSignalPOIs r_ggH --freezeParameters r_bbH,r_ggH \
        --X-rtd MINIMIZER_analytic \
        --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
        --robustHesse 1 \
        -n .combined-cmb.for_shape_unblinding \
        --there \
        -v 1

elif [[ $MODE == "postfit-plots" ]]; then
    ######################
    # Extract postfit shapes.
    #####################
    fitfile=${datacarddir}/combined/cmb/fitDiagnostics.combined-cmb.for_shape_unblinding.root
    for era in 2016 2017 2018; do
        prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/${era}/htt_em_2_*/combined.txt.cmb" \
                                          --workspace_name ws.root \
                                          --fit_arguments "-f ${fitfile}:fit_b --postfit --sampling" \
                                          --output_name postfit_shapes_${freeze}.root \
                                          --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log
        prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/${era}/htt_*_3*_*/combined.txt.cmb" \
                                          --workspace_name ws.root \
                                          --freeze_arguments "--freeze ${freeze}" \
                                          --fit_arguments "-f ${fitfile}:fit_b --postfit --sampling" \
                                          --output_name postfit_shapes_${freeze}.root \
                                          --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log
    done

    hadd -f ${datacarddir}/combined/cmb/postfit_shapes_${freeze}.root ${datacarddir}/201?/htt_*/postfit_shapes_${freeze}.root | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log

    for era in 2016 2017 2018; do
        bash plotting/plot_shapes_mssm_model_independent.sh \
            ${era} \
            "${datacarddir}/combined/cmb/postfit_shapes_${freeze}.root" \
            "${datacarddir}/plots/postfit_shapes_$(echo ${freeze} | sed 's/=//g; s/\./p/g')/" \
            et,mt,tt,em \
            $(echo $freeze | cut -d, -f1 | cut -d= -f2) \
            $(echo $freeze | cut -d, -f2 | cut -d= -f2)
    done

elif [[ $MODE == "prepare-ggH-bbH-scan" ]]; then
    [[ ! -d ${defaultdir}/ggH_bbH_scan_ind/condor ]] && mkdir -p ${defaultdir}/ggH_bbH_scan_ind/condor
    cd ${defaultdir}/ggH_bbH_scan_ind/condor
    # Run 2D likelihood scans for r_ggH and r_bbH
    combineTool.py -M MultiDimFit \
        --algo grid --points 225 --split-points 50 \
        -m "60,80,100,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500" \
        --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_ggH_bbH_2D_boundaries.json \
        --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_ggH,r_bbH \
        -d ${datacarddir}/combined/cmb/ws.root \
        --job-mode condor --dry-run --task-name ggH_bbH_likelihood_scan \
        --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
        -n ".ggH-bbH" \
        -v 1
        # --there -n ".ggH-bbH"
        # -m "60,80,100,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500" \

    # Create asimov dataset for SM-expectation.
    combineTool.py -M MultiDimFit \
        --algo none \
        -m "125" \
        --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_ggH_bbH_2D_boundaries.json \
        --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_ggH,r_bbH \
        -d ${datacarddir}/combined/cmb/ws.root \
        -t -1 --saveToys \
        --there -n ".2D.ToyDataset.SM1" \
        --dry-run
        # -d output/mssm_201017_SMHbkg/cmb/ws.root \ TODO: maybe need different datacard here

    # Copy it to correct location
    # Copy only necessary in case different datacards are used for the asimov creation and the fits.
    # cp output/mssm_070617_SMHbkg/cmb/higgsCombine.2D.ToyDataset.SM1.MultiDimFit.mH125.123456.root output/mssm_201017/cmb/higgsCombine.2D.ToyDataset.SM1.MultiDimFit.mH125.123456.root

    # Run fits on this asimov dataset
    combineTool.py -M MultiDimFit \
        --algo none \
        -m "60,80,100,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500" \
        --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_ggH_bbH_2D_boundaries.json \
        --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_ggH,r_bbH \
        -d ${datacarddir}/combined/cmb/ws.root \
        -t -1 --toysFile higgsCombine.2D.ToyDataset.SM1.MultiDimFit.mH125.123456.root \
        --job-mode condor --dry-run --task-name ggH_bbH_likelihood_SM --merge 3 \
        --there -n ".2D.SM1.bestfit"

elif [[ $MODE == "submit-ggH-bbH-scan" ]]; then
    cd ${defaultdir}/ggH_bbH_scan_ind/condor
    condor_submit condor_ggH_bbH_likelihood_scan.sub
    # condor_submit condor_ggH_bbH_likelihood_SM.sub

elif [[ $MODE == "collect-ggH-bbH-scan" ]]; then
    cd ${defaultdir}/ggH_bbH_scan_ind/
    for mass in 60 80 100 120 125 130 140 160 180 200 250 300 350 400 450 500 600 700 800 900 1000 1200 1400 1600 1800 2000 2300 2600 2900 3200 3500; do
        python ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/plotting/plotMultiDimFit.py \
            --title-right="138 fb^{-1} (13 TeV)" \
            --cms-sub="Preliminary" \
            --mass $mass \
            -o 2D_limit_mH${mass} \
            --debug-output test_mH${mass}.root \
            condor/higgsCombine.ggH-bbH.POINTS.*.mH${mass}.root
    done

fi
