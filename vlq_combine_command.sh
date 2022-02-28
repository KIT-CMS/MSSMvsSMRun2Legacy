ulimit -s unlimited

TAG=$1
MODE=$2
if [[ $TAG == "auto" ]]; then
    TAG="cmb_ind"
fi

analysis="vector_leptoquark"
hSM_treatment="hSM-in-bg"
categorization="classic"
sm_like_hists="bsm"
#sub_analyses="betaRd33_0 betaRd33_minus1 betaRd33_0_offdiag0"
sub_analyses="betaRd33_0"
#sub_analyses="betaRd33_minus1 betaRd33_0_offdiag0"
#sub_analyses="betaRd33_0_offdiag0"
MASS=$3
#gU=1.21
gU=2.52
group="--group htt_em --group htt_et --group htt_mt --group htt_tt"

for sub_analysis in $sub_analyses; do

  defaultdir=analysis/${TAG}_${sub_analysis}
  [[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
  [[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
  datacarddir=${defaultdir}/datacards_${analysis}


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
          --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_nobtag_categories.txt \
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
      -m ${MASS} --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

  elif [[ $MODE == "ws-no-int" ]]; then
      ############
      # workspace creation
      ############

      combineTool.py -M T2W -o "ws.root" \
      -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
      --PO no_interference \
      -i ${datacarddir}/combined/cmb/ \
      -m ${MASS} --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

  elif [[ $MODE == "ws-grid" ]]; then
      ############
      # workspace creation
      ############

      combineTool.py -M T2W -o "ws.root" \
      -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
      --PO grid \
      -i ${datacarddir}/combined/cmb/ \
      -m 1000 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

  elif [[ $MODE == "ws-ccc" ]]; then
      ############
      # workspace creation
      ############

      combineTool.py -M T2W -o "ws.root" \
      -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
      --PO mu \
      -i ${datacarddir}/combined/cmb/ \
      -m 1000 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

  elif [[ $MODE == "make-asimov" ]]; then
      combineTool.py -m "500,1000,2000,3000,4000,5000" \
      -M GenerateOnly -t -1 \
      --saveToys \
      --setParameters gU=0,lumi_scale=1,r_ggH=3 \
      -d ${datacarddir}/combined/cmb/ws.root \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      -n .$1 \
      --task-name vlq_${sub_analysis}_asimov_cmb 

 
  elif [[ $MODE == "submit" ]]; then 
      #combineTool.py -m "500,1000,2000,3000,4000,5000" \
      combineTool.py -m "2000" \
      -M AsymptoticLimits \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU \
      -d ${datacarddir}/combined/cmb/ws.root \
      --there -n ".limit" \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_full_cmb \
      --X-rtd MINIMIZER_analytic \
      --rAbsAcc 0 \
      --rRelAcc 0.0005 \
      --cminDefaultMinimizerStrategy 0 \
      --cminDefaultMinimizerTolerance 0.01 \
      -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh.txt

  elif [[ $MODE == "submit-grid" ]]; then
      cd ${datacarddir}/combined/cmb

      combineTool.py \
      -M AsymptoticGrid \
      ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/vlq_asymptotic_grid_test.json \
      --redefineSignalPOIs r \
      --setParameterRanges r=0,1 \
      --setParameters r=1,lumi_scale=1 \
      -d ws.root \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_full_cmb \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 \
      -v 1

  elif [[ $MODE == "collect-grid" ]]; then
      cd ${datacarddir}/combined/cmb

      combineTool.py \
      -M AsymptoticGrid \
      ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/vlq_asymptotic_grid_test.json \
      --redefineSignalPOIs r \
      --setParameterRanges r=0,1 \
      --setParameters r=1,lumi_scale=1 \
      -d ws.root \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_full_cmb_collect \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 \
      -v 1

      mv asymptotic_grid.root asymptotic_grid_${sub_analysis}.root

      cd ../../../../../
  elif [[ $MODE == "plot-grid" ]]; then
      if [[ ${sub_analysis} == "betaRd33_0" ]]; then
        scenario_label="VLQ benchmark 1"
      elif [[ ${sub_analysis} == "betaRd33_minus1" ]]; then
        scenario_label="VLQ benchmark 2"
      elif [[ ${sub_analysis} == "betaRd33_0_offdiag0" ]]; then
        scenario_label="VLQ benchmark 3"
      fi

      if [[ ${sub_analysis} != "betaRd33_0_offdiag0" ]]; then
        plotLimitGrid.py ${datacarddir}/combined/cmb/asymptotic_grid_${sub_analysis}.root \
        --title-left="${scenario_label}" \
        --cms-sub="" \
        --output vlq_${sub_analysis}_grid \
        --title-right="138 fb^{-1} (13 TeV)" \
        --contours="exp0,exp-1,exp-2,exp+1,exp+2,obs" \
        --y-range 0.0,6.0 \
        --x-title "m_{U} (TeV)" \
        --y-title "g_{U}" \
        --bin-method "BinCenterAligned" \
        --no_morph \
        --convert_GeV_to_TeV \
        --allowed_region_from_json "input/vlq_${sub_analysis}_bestfit.json" \
        --remove_contours_below "1"
      else
        plotLimitGrid.py ${datacarddir}/combined/cmb/asymptotic_grid_${sub_analysis}.root \
        --title-left="${scenario_label}" \
        --cms-sub="" \
        --output vlq_${sub_analysis}_grid \
        --title-right="138 fb^{-1} (13 TeV)" \
        --contours="exp0,exp-1,exp-2,exp+1,exp+2,obs" \
        --y-range 0.0,6.0 \
        --x-title "m_{U} (TeV)" \
        --y-title "g_{U}" \
        --bin-method "BinCenterAligned" \
        --no_morph \
        --convert_GeV_to_TeV \
        --remove_contours_below "1"

        #plotLimitGrid.py ${datacarddir}/combined/cmb/asymptotic_grid_${sub_analysis}.root \
        #--title-left="${scenario_label}" \
        #--cms-sub="" \
        #--output vlq_${sub_analysis}_grid \
        #--title-right="138 fb^{-1} (13 TeV)" \
        #--contours="exp0,exp-1,exp-2,exp+1,exp+2,obs" \
        #--y-range 0.0,6.0 \
        #--x-title "m_{U} (TeV)" \
        #--y-title "g_{U}" \
        #--bin-method "BinCenterAligned" \
        #--convert_GeV_to_TeV 
      fi



  elif [[ $MODE == "submit-extrap-run3" ]]; then
      # need a new directory for this
      combineTool.py -m "500,1000,2000,3000,4000,5000" \
      -M AsymptoticLimits \
      -t -1 \
      --setParameters gU=0,lumi_scale=3.62 \
      --redefineSignalPOIs gU \
      -d ${datacarddir}/combined/cmb/ws.root \
      --there -n ".limit" \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_full_cmb \
      --X-rtd MINIMIZER_analytic \
      --rAbsAcc 0 \
      --rRelAcc 0.0005 \
      --cminDefaultMinimizerStrategy 0 \
      --cminDefaultMinimizerTolerance 0.01 \
      -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh.txt

 
  elif [[ $MODE == "submit-extrap-hllhc" ]]; then
      # need a new directory for this
      combineTool.py -m "500,1000,2000,3000,4000,5000" \
      -M AsymptoticLimits \
      -t -1 \
      --setParameters gU=0,lumi_scale=21.74 \
      --redefineSignalPOIs gU \
      -d ${datacarddir}/combined/cmb/ws.root \
      --there -n ".limit" \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_full_cmb \
      --X-rtd MINIMIZER_analytic \
      --rAbsAcc 0 \
      --rRelAcc 0.0005 \
      --cminDefaultMinimizerStrategy 0 \
      --cminDefaultMinimizerTolerance 0.01 \
      -v 1 | tee -a ${defaultdir}/logs/job_setup_modelind_bbh.txt


  elif [[ $MODE == "significance" ]]; then
      ############
      # job setup creation
      ############
      combineTool.py -m "500,1000,2000,3000,4000,5000" \
      -M Significance \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU \
      -d ${datacarddir}/combined/cmb/ws.root \
      --there -n ".significance" \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_significance \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 \
      --cminDefaultMinimizerTolerance 0.01 \
      -v 1

      combineTool.py -m "500,1000,2000,3000,4000,5000" \
      -M Significance \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU \
      -d ${datacarddir}/combined/cmb/ws.root \
      --there -n ".pvalue" \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_pvalue \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 \
      --cminDefaultMinimizerTolerance 0.01 \
      --pvalue \
      -v 1

  elif [[ $MODE == "collect-only" ]]; then
      combineTool.py -M CollectLimits ${datacarddir}/combined/cmb/higgsCombine.limit*.root \
      --use-dirs \
      -o ${datacarddir}/combined/cmb/vlq_${sub_analysis}.json
 
  
  elif [[ $MODE == "collect" ]]; then
      combineTool.py -M CollectLimits ${datacarddir}/combined/cmb/higgsCombine.limit*.root \
      --use-dirs \
      -o ${datacarddir}/combined/cmb/vlq_${sub_analysis}.json
   
      #plotMSSMLimits_backup.py --cms-sub "Preliminary" \
      #--title-right "138 fb^{-1} (13 TeV)" \
      #--x-title "M_{U} [TeV]"\
      #--y-axis-min 0.0 \
      #--y-axis-max 6.0 \
      #--show obs,exp ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
      #--process "vector_leptoquark" \
      #--subprocess "${sub_analysis}" \
      #--add-exp-line-from-json "{\"EXO-19-016 Total Expected\":\"input/total_expected_kappa_1.json\",\"EXO-19-016 Nonres Expected\":\"input/nonres_expected_kappa_1.json\"}" \
      #--add-obs-line-from-json "{\"EXO-19-016 Total Observed\":\"input/total_observed_kappa_1.json\",\"EXO-19-016 Nonres Observed\":\"input/nonres_observed_kappa_1.json\"}" \
      #--output vlq_${sub_analysis}_cmb_exo_kapp_1

      #plotMSSMLimits.py --cms-sub "Preliminary" \
      #--title-right "138 fb^{-1} (13 TeV)" \
      #--x-title "M_{U} [TeV]"\
      #--y-axis-min 0.0 \
      #--y-axis-max 6.0 \
      #--show obs,exp ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
      #--process "vector_leptoquark" \
      #--subprocess "${sub_analysis}" \
      #--add-exp-line-from-json "{\"3000 fb^{-1}\":\"analysis/1911_HLLHC_no_interference_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\",\"500 fb^{-1}\":\"analysis/1911_Run3_no_interference_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\"}" \
      #--output vlq_${sub_analysis}_cmb


      #plotMSSMLimits.py --cms-sub "Preliminary" \
      #--title-right "138 fb^{-1} (13 TeV)" \
      #--x-title "M_{U} [TeV]"\
      #--y-axis-min 0.0 \
      #--y-axis-max 6.0 \
      #--show obs,exp ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
      #--process "vector_leptoquark" \
      #--subprocess "${sub_analysis}" \
      #--add-exp-line-from-json "{\"3000 fb^{-1}\":\"analysis/2610_HLLHC_lumiscale_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\",\"500 fb^{-1}\":\"analysis/2610_run3_lumiscale_v2_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\"}" \
      #--output vlq_${sub_analysis}_cmb

      plotMSSMLimits_backup.py --cms-sub "Preliminary" \
      --title-right "138 fb^{-1} (13 TeV)" \
      --x-title "M_{U} [TeV]"\
      --y-axis-min 0.0 \
      --y-axis-max 6.0 \
      --show obs,exp0 ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
      --process "vector_leptoquark" \
      --subprocess "${sub_analysis}" \
      --add-exp-line-from-json "{\"3000 fb^{-1}\":\"analysis/1901_extrap_hllhc_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\",\"500 fb^{-1}\":\"analysis/1901_extrap_run3_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\"}" \
      --output vlq_${sub_analysis}_cmb_extrap


      plotMSSMLimits_backup.py --cms-sub "Preliminary" \
      --title-right "138 fb^{-1} (13 TeV)" \
      --x-title "M_{U} [TeV]"\
      --y-axis-min 0.0 \
      --y-axis-max 6.0 \
      --show obs,exp ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
      --process "vector_leptoquark" \
      --subprocess "${sub_analysis}" \
      --output vlq_${sub_analysis}_cmb_no_extrap


#
#      plotMSSMLimits.py --cms-sub "Preliminary" \
#      --title-right "138 fb^{-1} (13 TeV)" \
#      --x-title "M_{U} [TeV]"\
#      --y-axis-min 0.0 \
#      --y-axis-max 4.5 \
#      --show exp ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
#      --process "vector_leptoquark" \
#      --subprocess "${sub_analysis}" \
#      --add-exp-line-from-json "{\"500 fb^{-1}\":\"analysis/2610_run3_lumiscale_v2_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\",\"3000 fb^{-1}\":\"analysis/2610_HLLHC_lumiscale_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\"}" \
#      --convert-gU-to_lambda \
#      --output vlq_${sub_analysis}_cmb_lambda


  elif [[ $MODE == "ws-grid" ]]; then
      ############
      # workspace creation
      ############

      combineTool.py -M T2W -o "ws.root" \
      --PO grid \
      -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
      -i ${datacarddir}/combined/cmb/ \
      -m 1000 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

  elif [[ $MODE == "run-impacts" ]]; then
      taskname="impacts_${TAG}_${sub_analysis}_mH${MASS}"

      combineTool.py -M Impacts -d ${datacarddir}/combined/cmb/ws.root \
      --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
      --doInitialFit --robustFit 1 \
      -m $MASS \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU  -v 0

      #combineTool.py -M Impacts -d ${datacarddir}/combined/cmb/ws.root \
      #--X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
      #--robustFit 1 --doFits \
      #-m $MASS \
      #--setParameters gU=0,lumi_scale=1 \
      #--redefineSignalPOIs gU \
      #--job-mode "SGE" \
      #--prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      #--task-name ${taskname} --merge 5


  elif [[ $MODE == "signal-strength" ]]; then

      ulimit -s unlimited

      combineTool.py -M FitDiagnostics \
      -d ${datacarddir}/combined/cmb/ws.root \
      -m $MASS \
      --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU  -v 0 \
      --robustFit 1 --setRobustFitAlgo Minuit2 --setRobustFitStrategy 0 \ 
      --setRobustFitTolerance 0.2 \ 
      --cminPreScan --cminPreFit 1 \
      --rMin=0 --rMax=4


  elif [[ $MODE == "collect-impacts" ]]; then
      combineTool.py -M Impacts -d ${datacarddir}/combined/cmb/ws.root -m $MASS -o ${sub_analysis}_${MASS}_impacts.json --redefineSignalPOIs gU

  elif [[ $MODE == "plot-impacts" ]]; then
      plotImpacts.py -i ${sub_analysis}_${MASS}_impacts.json -o ${sub_analysis}_${MASS}_impacts --transparent

  elif [[ $MODE == "ws-plot" ]]; then

      combineTool.py -M T2W -o "ws.root" \
      -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
      -i ${datacarddir}/201?/htt_*/ \
      --X-allow-no-signal \
      -m ${MASS} --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt


  elif [[ $MODE == "fit-for-plots" ]]; then
      combineTool.py -M FitDiagnostics \
      -d ${datacarddir}/combined/cmb/ws.root \
      -m $MASS \
      --setParameters gU=$gU,lumi_scale=1 \
      --redefineSignalPOIs gU --freezeParameters gU \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
      --robustHesse 1 \
      -n .combined-cmb.for_shape_unblinding \
      --there \
      -v 1

  elif [[ $MODE == "bkg-only-fit-for-plots" ]]; then
      combineTool.py -M FitDiagnostics \
      -d ${datacarddir}/combined/cmb/ws.root \
      -m $MASS \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU --freezeParameters gU \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
      --robustHesse 1 \
      -n .combined-cmb.for_shape_unblinding_bkg \
      --there \
      -v 1

  elif [[ $MODE == "fit-for-scan" ]]; then

      cd ${datacarddir}/combined/cmb

      #combine -M MultiDimFit \
      #-m $MASS \
      #--setParameters gU=0.8885,lumi_scale=3.62  \
      #-t -1 --freezeParameters gU \
      #--redefineSignalPOIs gU \
      #--X-rtd MINIMIZER_analytic \
      #--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
      #--robustHesse 1 \
      #-d ${datacarddir}/combined/cmb/ws.root \
      #--points=30 --fastScan --algo=grid --rMin=0 --rMax=3 \
      #-n "Scan"

      combine -M MultiDimFit \
      -m $MASS \
      --redefineSignalPOIs gU \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
      --robustHesse 1 \
      -d ws.root \
      --points=30 --fastScan --algo=grid --rMin=0 --rMax=2 \
      -n "Scan"

      cd ../../../../../

  elif [[ $MODE == "fit-for-scan-asimov" ]]; then

      cd ${datacarddir}/combined/cmb

      #combine -M MultiDimFit \
      #-m $MASS \
      #--setParameters gU=0.8885,lumi_scale=3.62  \
      #-t -1 --freezeParameters gU \
      #--redefineSignalPOIs gU \
      #--X-rtd MINIMIZER_analytic \
      #--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
      #--robustHesse 1 \
      #-d ${datacarddir}/combined/cmb/ws.root \
      #--points=30 --fastScan --algo=grid --rMin=0 --rMax=3 \
      #-n "Scan"

      if [[ $MASS == "500" ]]; then
        combine -M MultiDimFit \
        -m $MASS \
        --toysFile="../../../../../higgsCombine.${1}.GenerateOnly.mH${MASS}.123456.root" \
        --redefineSignalPOIs gU \
        --setParameters r_ggH=0  \
        --freezeParameters r_ggH \
        --X-rtd MINIMIZER_analytic \
        --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
        --robustHesse 1 \
        -d ws.root \
        --points=30 --fastScan --algo=grid --rMin=0 --rMax=2 \
        -n "Scan"
      elif [[ $MASS == "2000" ]]; then
        combine -M MultiDimFit \
        -m $MASS \
        --toysFile="../../../../../higgsCombine.${1}.GenerateOnly.mH${MASS}.123456.root" \
        --redefineSignalPOIs gU \
        --setParameters r_ggH=0  \
        --freezeParameters r_ggH \
        --X-rtd MINIMIZER_analytic \
        --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
        --robustHesse 1 \
        -d ws.root \
        --points=30 --fastScan --algo=grid --rMin=0 --rMax=4 \
        -n "Scan"
      fi


      cd ../../../../../


  elif [[ $MODE == "lkld-plot" ]]; then
    ../CombineTools/scripts/plot1DScan.py ${datacarddir}/combined/cmb/higgsCombineScan.MultiDimFit.mH${MASS}.root --POI gU  --y-cut=8 --output=${datacarddir}/combined/cmb/scan_$MASS


  elif [[ $MODE == "prefit-plots" ]]; then
      freeze="MH=${MASS},gU=${gU}"
      for era in 2016 2017 2018; do
          #prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/${era}/htt_em_2_*/combined.txt.cmb" \
          #                                  --workspace_name ws.root \
          #                                  --output_name prefit_shapes_${freeze}.root \
          #                                  --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-combined-${freeze}.log
          prefit_postfit_shapes_parallel.py --datacard_pattern "${datacarddir}/${era}/htt_*_3*/combined.txt.cmb" \
                                            --workspace_name ws.root \
                                            --freeze_arguments "--freeze ${freeze}" \
                                            --output_name prefit_shapes_${freeze}.root \
                                            --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-combined-${freeze}.log
      done
      hadd -f ${datacarddir}/combined/cmb/prefit_shapes_${freeze}.root ${datacarddir}/201?/htt_*/prefit_shapes_${freeze}.root | tee -a ${defaultdir}/logs/extract_model_independent_shapes-combined-${freeze}.log

      for era in 2016 2017 2018; do
          bash plotting/plot_shapes_vlq.sh \
              ${era} \
              "${datacarddir}/combined/cmb/prefit_shapes_${freeze}.root" \
              "${datacarddir}/plots/prefit_shapes_$(echo ${freeze} | sed 's/=//g; s/\./p/g')/" \
              tt,mt,et,em \
              $MASS \
              $gU \
              $sub_analysis
      done

  elif [[ $MODE == "postfit-plots" ]]; then
      freeze="MH=${MASS},gU=${gU}"
      fitfile=${datacarddir}/combined/cmb/fitDiagnostics.combined-cmb.for_shape_unblinding.root
      echo $fitfile
      for era in 2016 2017 2018; do
          eval "prefit_postfit_shapes_parallel.py --datacard_pattern \"${datacarddir}/${era}/htt_em_2_*/combined.txt.cmb\" \
                                            --workspace_name ws.root \
                                            --fit_arguments \"-f ${fitfile}:fit_b --postfit --sampling\" \
                                            --output_name postfit_shapes_${freeze}.root \
                                            --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log"
          eval "prefit_postfit_shapes_parallel.py --datacard_pattern \"${datacarddir}/${era}/htt_*_3*_*/combined.txt.cmb\" \
                                            --workspace_name ws.root \
                                            --freeze_arguments \"--freeze ${freeze}\" \
                                            --fit_arguments \"-f ${fitfile}:fit_b --postfit --sampling\" \
                                            --output_name postfit_shapes_${freeze}.root \
                                            --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log"
      done
  
      hadd -f ${datacarddir}/combined/cmb/postfit_shapes_${freeze}.root ${datacarddir}/201?/htt_*/postfit_shapes_${freeze}.root | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log
 
      mkdir ${datacarddir}/plots 
      for era in 2016 2017 2018; do
          bash plotting/plot_shapes_vlq.sh \
              ${era} \
              "${datacarddir}/combined/cmb/postfit_shapes_${freeze}.root" \
              "${datacarddir}/plots/postfit_shapes_$(echo ${freeze} | sed 's/=//g; s/\./p/g')/" \
              et,mt,tt,em \
              $MASS \
              $gU \
              $sub_analysis
      done

  elif [[ $MODE == "postfit-plots-bkg" ]]; then
      freeze="MH=${MASS},gU=${gU}"
      fitfile=${datacarddir}/combined/cmb/fitDiagnostics.combined-cmb.for_shape_unblinding_bkg.root
      echo $fitfile
      for era in 2016 2017 2018; do
          eval "prefit_postfit_shapes_parallel.py --datacard_pattern \"${datacarddir}/${era}/htt_em_2_*/combined.txt.cmb\" \
                                            --workspace_name ws.root \
                                            --fit_arguments \"-f ${fitfile}:fit_b --postfit --sampling\" \
                                            --output_name postfit_shapes_bkg_${freeze}.root \
                                            --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log"
          eval "prefit_postfit_shapes_parallel.py --datacard_pattern \"${datacarddir}/${era}/htt_*_3*_*/combined.txt.cmb\" \
                                            --workspace_name ws.root \
                                            --freeze_arguments \"--freeze ${freeze}\" \
                                            --fit_arguments \"-f ${fitfile}:fit_b --postfit --sampling\" \
                                            --output_name postfit_shapes_bkg_${freeze}.root \
                                            --parallel 8 | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log"
      done

      hadd -f ${datacarddir}/combined/cmb/postfit_shapes_bkg_${freeze}.root ${datacarddir}/201?/htt_*/postfit_shapes_bkg_${freeze}.root | tee -a ${defaultdir}/logs/extract_model_independent_shapes-postfit-combined-${freeze}.log

      mkdir ${datacarddir}/plots
      for era in 2016 2017 2018; do
          bash plotting/plot_shapes_vlq.sh \
              ${era} \
              "${datacarddir}/combined/cmb/postfit_shapes_bkg_${freeze}.root" \
              "${datacarddir}/plots/postfit_shapes_bkg_$(echo ${freeze} | sed 's/=//g; s/\./p/g')/" \
              et,mt,tt,em \
              $MASS \
              $gU \
              $sub_analysis
      done


  elif [[ $MODE == "ccc-channel" ]]; then

      ulimit -s unlimited

      combine -M ChannelCompatibilityCheckRegexGroup \
      -d ${datacarddir}/combined/cmb/ws.root \
      -m $MASS \
      --setParameters lumi_scale=1 \
      --setParameterRanges mu=-100,30 \
      --redefineSignalPOIs mu \
      --X-rtd MINIMIZER_analytic \
      --X-rtd FITTER_NEW_CROSSING_ALGO \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
      --robustFit 1 --stepSize 0.01 \
      -n .CCC.channel \
      --saveFitResult \
      -g em -g et -g mt -g tt \
      -v 1

  elif [[ $MODE == "ccc-year" ]]; then

      ulimit -s unlimited

      combine -M ChannelCompatibilityCheckRegexGroup \
      -d ${datacarddir}/combined/cmb/ws.root \
      -m $MASS \
      --setParameters lumi_scale=1 \
      --setParameterRanges mu=-2,6 \
      --redefineSignalPOIs mu \
      --X-rtd MINIMIZER_analytic \
      --X-rtd FITTER_NEW_CROSSING_ALGO \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
      --robustFit 1 --stepSize 0.01 \
      -n .CCC.year \
      --saveFitResult \
      -g 2016 -g 2017 -g 2018 \
      -v 1

  elif [[ $MODE == "ccc-btag" ]]; then

      ulimit -s unlimited

      combine -M ChannelCompatibilityCheckRegexGroup \
      -d ${datacarddir}/combined/cmb/ws.root \
      -m $MASS \
      --setParameters lumi_scale=1 \
      --setParameterRanges mu=-2,20 \
      --redefineSignalPOIs mu \
      --X-rtd MINIMIZER_analytic \
      --X-rtd FITTER_NEW_CROSSING_ALGO \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
      --robustFit 1 --stepSize 0.01 \
      -n .CCC.btag \
      --saveFitResult \
      -g btag -g nobtag  \
      -v 1



  elif [[ $MODE == "plot-ccc-channel" ]]; then
      python plotting/plot_ccc.py \
      higgsCombine.CCC.channel.ChannelCompatibilityCheckRegexGroup.mH${MASS}.root \
      -o ChannelCompatibilityCheck_FitResults_mH${MASS}_channel \
      -p mu \
      -r m60,60

  elif [[ $MODE == "plot-ccc-year" ]]; then
      python plotting/plot_ccc.py \
      higgsCombine.CCC.year.ChannelCompatibilityCheckRegexGroup.mH${MASS}.root \
      -o ChannelCompatibilityCheck_FitResults_mH${MASS}_year \
      -p mu \
      -r m6,14

  elif [[ $MODE == "plot-ccc-btag" ]]; then
      python plotting/plot_ccc.py \
      higgsCombine.CCC.btag.ChannelCompatibilityCheckRegexGroup.mH${MASS}.root \
      -o ChannelCompatibilityCheck_FitResults_mH${MASS}_btag \
      -p mu \
      -r m10,35

  fi
done
