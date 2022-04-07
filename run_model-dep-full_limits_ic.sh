# based on https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/ntuple_processor
ulimit -s unlimited

TAG=$1
MODE=$2
MODEL=$3
ANALYSISTYPE=$4
HSMTREATMENT=$5
OLDFILES=$6
[[ -z $6 ]] && OLDFILES=0

if [[ $ANALYSISTYPE == "classic" ]]; then
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="classic"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_classic"
    fi
elif [[ $ANALYSISTYPE == "lowmass" ]]; then
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="lowmass"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_lowmass"
    fi
elif [[ $ANALYSISTYPE == "couplings" ]]; then
    analysis="bsm-model-dep-additional"
    sm_like_hists="sm125"
    categorization="lowmass"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_couplings"
    fi
else
    analysis="bsm-model-dep-full"
    sm_like_hists="sm125"
    replace_with_sm125=1
    categorization="with-sm-ml"
    if [[ $TAG == "auto" ]]; then
        TAG="cmb_with_ml"
    fi
fi
# Szenarios from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWGMSSMNeutral?redirectedfrom=LHCPhysics.LHCHXSWGMSSMNeutral#Baseline_scenarios
### MSSM scenarios #####
scale_qqh_by_hand=0
if [[ $MODEL == "mh125" ]]; then
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    scenario_label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_lc" ]]; then
    wsoutput="ws_mh125_lc.root"
    modelfile="13,Run2017,mh125_lc_13.root"
    scenario_label="M_{h}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_ls" ]]; then
    wsoutput="ws_mh125_ls.root"
    modelfile="13,Run2017,mh125_ls_13.root"
    scenario_label="M_{h}^{125}(#tilde{#tau}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
elif [[ $MODEL == "mh125_align" ]]; then
    wsoutput="ws_mh125_align.root"
    modelfile="13,Run2017,mh125_align_13.root"
    scenario_label="M_{h}^{125}(alignment) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=3.0
    y_max=12.0
elif [[ $MODEL == "mHH125" ]]; then
    wsoutput="ws_mHH125.root"
    modelfile="13,Run2017,mHH125_13.root"
    scenario_label="M_{H}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-heavy"
    sm_like_mass="m_H"
    x_title='m_{H^{#plus}} (GeV)'
    mass_histogram_title="m_{H}"
    y_min=5.0
    y_max=6.0
elif [[ $MODEL == "mh1125_CPV" ]]; then
    wsoutput="ws_mh1125_cpv.root"
    modelfile="13,Run2017,mh1125_CPV_13.root"
    scenario_label="M_{h_{1}}^{125}(CPV) scenario (^{}h_{1},^{}h_{2},^{}h_{3}#rightarrow#tau#tau)"
    sub_analysis="cpv"
    sm_like_mass="m_H1"
    x_title='m_{H^{#plus}} (GeV)'
    mass_histogram_title="m_{^{}h_{1}}"
    y_min=1.0
    y_max=20.0
### Negative mu scenarios #####
elif [[ $MODEL == "mh125_muneg_1" ]]; then
    wsoutput="mh125_muneg_1.root"
    modelfile="13,Run2017,mh125_muneg_1_13.root"
    scenario_label="M_{h}^{125 ^{}#mu_{1}#minus} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=4.0
    y_max=56.0
elif [[ $MODEL == "mh125_muneg_2" ]]; then
    wsoutput="mh125_muneg_2.root"
    modelfile="13,Run2017,mh125_muneg_2_13.root"
    scenario_label="M_{h}^{125 ^{}#mu_{2}#minus} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=5.0
    y_max=30.0
elif [[ $MODEL == "mh125_muneg_3" ]]; then
    wsoutput="mh125_muneg_3.root"
    modelfile="13,Run2017,mh125_muneg_3_13.root"
    scenario_label="M_{h}^{125 ^{}#mu_{3}#minus} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=6.0
    y_max=20.0
### EFT scenarios #####
elif [[ $MODEL == "mh125EFT" ]]; then
    wsoutput="mh125EFT.root"
    modelfile="13,Run2017,mh125EFT_13.root"
    scenario_label="M_{h,EFT}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=10.0
elif [[ $MODEL == "mh125EFT_lc" ]]; then
    wsoutput="mh125EFT_lc.root"
    modelfile="13,Run2017,mh125EFT_lc_13.root"
    scenario_label="M_{h,EFT}^{125}(#tilde{#chi}) scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=10.0
### hMSSM scenario #####
elif [[ $MODEL == "hMSSM" ]]; then
    wsoutput="hMSSM.root"
    modelfile="13,Run2017,hMSSM_13.root"
    scenario_label="hMSSM scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
else
    wsoutput="ws_mh125.root"
    modelfile="13,Run2017,mh125_13.root"
    scenario_label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)"
    sub_analysis="sm-like-light"
    sm_like_mass="m_h"
    x_title='m_{A} (GeV)'
    mass_histogram_title="m_{h}"
    y_min=1.0
    y_max=60.0
fi

if [[ $MODEL == "couplings" ]]; then
    wsoutput="ws_couplings.root"
    sub_analysis="sm-like-light"
fi

if [[ $OLDFILES == 1 ]]; then
    wsoutput=${wsoutput/.root/_old.root}
    modelfile=${modelfile/.root/_old.root}
fi

defaultdir="analysis/$TAG"
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
defaultdir=$(readlink -f ${defaultdir})
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/limits_${MODEL}/jobs ]] && mkdir -p ${defaultdir}/limits_${MODEL}/jobs

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
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log.txt
    elif [[ $ANALYSISTYPE == "lowmass" || $ANALYSISTYPE == "couplings" ]]; then

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 " \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_2d_to_1d.txt \
            --variable m_sv_VS_pt_tt_splitpT \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log_nobtag.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 " \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_1d_btag.txt \
            --variable m_sv_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log_btag.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 " \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_cr.txt \
            --variable mt_tot_puppi \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_mssm_log_ttbar.txt


    else
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_neuralnet_categories.txt \
            --variable nnscore \
            --sm \
            --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_sm_log.txt

        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis ${sub_analysis} \
            --hSM-treatment $HSMTREATMENT  \
            --categorization ${categorization} \
            --sm-like-hists ${sm_like_hists} \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
            --additional-arguments "--auto_rebin=1 --real_data=1 --manual_rebin=1 --split_sm_signal_cat=1 --enable_bsm_lowmass=1" \
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
    # Check if the expected number of datacards has been written
    EXPECTED=$(((15+15+11+18)*3))
    if [[ $(ls ${datacarddir}/combined/cmb/*.txt | wc -l) != $EXPECTED ]]; then
        echo "[ERROR] Not all datacards have been created or written. Please check the logs..."
        echo "Expected ${EXPECTED} datacards written but found only $(ls ${datacarddir}/combined/cmb/ | wc -l) in the combined directory."
    fi

elif [[ $MODE == "ws" ]]; then
    ############
    # workspace creation
    ############
    if [[ $ANALYSISTYPE == "couplings" ]]; then
        combineTool.py -M T2W -o ${wsoutput} \
        -P CombineHarvester.MSSMvsSMRun2Legacy.YtYbScan:YtYbScan \
        --PO XS-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/xs_lowmass_yb_yt.root \
        -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt
    elif [[ $OLDFILES == 0 ]]; then
        combineTool.py -M T2W -o ${wsoutput} \
        -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM \
        --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
        --PO replace-with-SM125=${replace_with_sm125} \
        --PO hSM-treatment=$HSMTREATMENT \
        --PO modelFile=${modelfile} \
        --PO minTemplateMass=60 \
        --PO maxTemplateMass=3500 \
        --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --PO sm-predictions=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_predictions_13TeV.json \
        --PO qqh-pred-from-scaling=${scale_qqh_by_hand} \
        -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt
    else
        combineTool.py -M T2W -o ${wsoutput} \
        -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM_oldModels:MSSMvsSM_oldModels \
        --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ \
        --PO replace-with-SM125=${replace_with_sm125} \
        --PO hSM-treatment=$HSMTREATMENT \
        --PO modelFile=${modelfile} \
        --PO minTemplateMass=60 \
        --PO maxTemplateMass=3500 \
        --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2.root \
        --PO sm-predictions=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_predictions_13TeV.json \
        -i ${datacarddir}/combined/cmb/ 2>&1 | tee -a ${defaultdir}/logs/workspace_${MODEL}.txt
    fi
elif [[ $MODE == "submit" && $ANALYSISTYPE == "couplings" ]]; then
  ############
  # job setup creation
  ############
  cd ${defaultdir}/limits_${MODEL}/jobs
  MASS="100"
  for MASS in "100"; do 

    combineTool.py -m "${MASS}" -M MultiDimFit --setParameterRanges Yt_H=0,0.7:Yb_H=-2.0,2.0 --freezeParameters r,Yt_A,Yb_A --setParameters mA=${MASS},mH=${MASS},r=1,Yt_H=0,Yb_H=0,Yt_A=0,Yb_A=0 --redefineSignalPOIs Yt_H,Yb_H \
  -d ${datacarddir}/combined/cmb/${wsoutput} \
  -n ".Yt_H_vs_Yb_H" --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
  --job-mode 'SGE' \
  --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
  --task-name "Ys_scan_mH${MASS}" \
  --algo grid --points 4000 --split-points 5 --alignEdges 1

    combineTool.py -m "${MASS}" -M MultiDimFit --setParameterRanges Yt_A=0,0.7:Yb_A=-2.0,2.0 --freezeParameters r,Yt_H,Yb_H --setParameters mA=${MASS},mH=${MASS},r=1,Yt_H=0,Yb_H=0,Yt_A=0,Yb_A=0 --redefineSignalPOIs Yt_A,Yb_A \
  -d ${datacarddir}/combined/cmb/${wsoutput} \
  -n ".Yt_A_vs_Yb_A" --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
  --job-mode 'SGE' \
  --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
  --task-name "Ys_scan_mA${MASS}" \
  --algo grid --points 4000 --split-points 5 --alignEdges 1

  done

elif [[ $MODE == "submit" ]]; then

    ############
    # job setup creation
    ############
    cd ${defaultdir}/limits_${MODEL}/jobs
    if [[ $HSMTREATMENT == "hSM-in-bg" ]]; then
        combineTool.py -M AsymptoticGrid \
        ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_${MODEL}.json \
        -d ${datacarddir}/combined/cmb/${wsoutput} \
        --job-mode 'SGE' \
        --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
        --task-name $taskname \
        --redefineSignalPOI r \
        --setParameterRanges r=0,1 \
        --setParameters r=1,x=1 \
        --freezeParameters x -v1 \
        --cminDefaultMinimizerStrategy 0 \
        --X-rtd MINIMIZER_analytic \
        --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}.txt
    elif [[ $HSMTREATMENT == "no-hSM-in-bg" ]]; then
        combineTool.py -M AsymptoticGrid \
        ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_${MODEL}.json \
        -d ${datacarddir}/combined/cmb/${wsoutput} \
        --job-mode 'SGE' \
        --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
        --task-name $taskname \
        --redefineSignalPOI x \
        --setParameterRanges x=0,1 \
        --setParameters r=1 \
        --freezeParameters r -v1 \
        --cminDefaultMinimizerStrategy 0 \
        --X-rtd MINIMIZER_analytic \
        --cminDefaultMinimizerTolerance 0.01 2>&1 | tee -a ${defaultdir}/logs/job_setup_${MODEL}.txt
    fi

elif [[ $MODE == "resubmit" ]]; then

   cd ${defaultdir}/limits_${MODEL}/jobs
   # if jog does exist then resubmit
   for job in $(ls *.sh); do
     log=${job/.sh/_%J.log}
     if test -f "${log}"; then continue; fi
     echo Log missing, resubmitting $job 
     #qsub -o ${defaultdir}/limits_${MODEL}/jobs/$log -q hep.q -l h_rt=3:0:0 $job
     qsub -o ${defaultdir}/limits_${MODEL}/jobs/$log -q hep.q -l h_rt=10:0:0 $job
   done
   # if log exists but it didn't finish then resubmit
   for out in $(grep -irL "Done in" *.log); do 
     job=${out/_%J.log/.sh}
     echo resubmitting $job 
     #qsub -o ${defaultdir}/limits_${MODEL}/jobs/$out -q hep.q -l h_rt=3:0:0 $job
     qsub -o ${defaultdir}/limits_${MODEL}/jobs/$out -q hep.q -l h_rt=10:0:0 $job
   done 

elif [[ $MODE == "collect" ]]; then
    ############
    # job collection
    ############
    cd ${defaultdir}/limits_${MODEL}/jobs
    combineTool.py -M AsymptoticGrid \
    ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_${MODEL}.json \
    -d ${datacarddir}/combined/cmb/${wsoutput} \
    --redefineSignalPOI x \
    --setParameterRanges x=0,1 \
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
    if [[ $ANALYSISTYPE == "classic" ]]; then
        title="Classic categorisation 138 fb^{-1} (13 TeV)"
    else
        title="138 fb^{-1} (13 TeV)"
    fi
    modelname=${MODEL}_13.root
    [[ $OLDFILES == 1 ]] && modelname="${MODEL}_13_old.root"
    for label in "Preliminary" ""; do
        plotLimitGrid.py asymptotic_grid.root \
        --scenario-label="${scenario_label}" \
        --output ${TAG}_${MODEL}_${label} \
        --title-right="${title}" \
        --cms-sub=${label} \
        --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
        --y-range ${y_min},${y_max} \
        --mass_histogram ${sm_like_mass} \
        --mass_histogram_title ${mass_histogram_title} \
        --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${modelname} \
        --x-title "${x_title}" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt
    done
elif [[ $MODE == "plot" ]]; then
    cd ${defaultdir}/limits_${MODEL}/

    ############
    # limit plot
    ############
    if [[ $ANALYSISTYPE == "classic" ]]; then
        title="Classic categorisation 138 fb^{-1} (13 TeV)"
    else
        title="138 fb^{-1} (13 TeV)"
    fi
    modelname=${MODEL}_13.root
    [[ $OLDFILES == 1 ]] && modelname="${MODEL}_13_old.root"
    for label in "Preliminary" ""; do
        ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/scripts/plotLimitGridMod.py asymptotic_grid.root \
        --scenario-label="${scenario_label}" \
        --output ${TAG}_${MODEL}_${label} \
        --title-right="${title}" \
        --cms-sub=${label} \
        --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" \
        --y-range ${y_min},${y_max} \
        --mass_histogram ${sm_like_mass} \
        --mass_histogram_title ${mass_histogram_title} \
        --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/${modelname} \
#        --extra_contour_file=/vols/cms/dw515/MSSMLowMass/CMSSW_10_2_25/src/CombineHarvester/MSSMvsSMRun2Legacy/MH125_contours.root \
#        --extra_contour_style=2,2,2,2,2,2,2 \
        --x-title "${x_title}" 2>&1 | tee -a ${defaultdir}/logs/plot_grid_${MODEL}.txt 
    done

fi
