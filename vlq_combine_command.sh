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
sub_analyses="betaRd33_0 betaRd33_minus1 betaRd33_0_offdiag0"
sub_analyses="betaRd33_minus1"
MASS=500
#gU=1.21
gU=0.8
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
          --eras 2018 \
          --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_tt_categories.txt \
          --additional-arguments "--auto_rebin=1 --real_data=1" \
          --variable mt_tot_puppi \
          --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_vlq_log.txt
 
      #morph_parallel.py --output ${defaultdir}/datacards \
      #    --analysis ${analysis} \
      #    --sub-analysis ${sub_analysis} \
      #    --hSM-treatment ${hSM_treatment} \
      #    --categorization ${categorization} \
      #    --sm-like-hists ${sm_like_hists} \
      #    --eras 2018 \
      #    --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_tt_categories.txt \
      #    --additional-arguments "--auto_rebin=0 --real_data=1" \
      #    --variable mt_tot_puppi \
      #    --parallel 10 2>&1 | tee -a ${defaultdir}/logs/morph_vlq_log.txt

 
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
      -m 1000 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

      #combineTool.py -M T2W -o "ws.root" \
      #-P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
      #--PO no_interference \
      #-i ${datacarddir}/combined/cmb/ \
      #-m 1000 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt

 
  elif [[ $MODE == "submit" ]]; then 
      ############
      # job setup creation
      ############
      combineTool.py -m "500,1000,2000,3000,4000,5000" \
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

 
  
  elif [[ $MODE == "collect" ]]; then
      combineTool.py -M CollectLimits ${datacarddir}/combined/cmb/higgsCombine.limit*.root \
      --use-dirs \
      -o ${datacarddir}/combined/cmb/vlq_${sub_analysis}.json
   
      #plotMSSMLimits.py --cms-sub "Preliminary" \
      #--title-right "138 fb^{-1} (13 TeV)" \
      #--x-title "M_{U} [TeV]"\
      #--y-axis-min 0.0 \
      #--y-axis-max 6.0 \
      #--show obs,exp ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
      #--process "vector_leptoquark" \
      #--subprocess "${sub_analysis}" \
      #--add-exp-line-from-json "{\"EXO-19-016 Total Expected\":\"input/total_expected_kappa_0.json\",\"EXO-19-016 Nonres Expected\":\"input/nonres_expected_kappa_0.json\"}" \
      #--add-obs-line-from-json "{\"EXO-19-016 Total Observed\":\"input/total_observed_kappa_0.json\",\"EXO-19-016 Nonres Observed\":\"input/nonres_observed_kappa_0.json\"}" \
      #--output vlq_${sub_analysis}_cmb

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

      plotMSSMLimits.py --cms-sub "Preliminary" \
      --title-right "138 fb^{-1} (13 TeV)" \
      --x-title "M_{U} [TeV]"\
      --y-axis-min 0.0 \
      --y-axis-max 6.0 \
      --show obs,exp0 ${datacarddir}/combined/cmb/vlq_${sub_analysis}_cmb.json \
      --process "vector_leptoquark" \
      --subprocess "${sub_analysis}" \
      --add-exp-line-from-json "{\"3000 fb^{-1}\":\"analysis/2610_HLLHC_lumiscale_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\",\"500 fb^{-1}\":\"analysis/2610_run3_lumiscale_v2_${sub_analysis}/datacards_vector_leptoquark/combined/cmb/vlq_${sub_analysis}_cmb.json\"}" \
      --output vlq_${sub_analysis}_cmb_extrap


      plotMSSMLimits.py --cms-sub "Preliminary" \
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

  elif [[ $MODE == "submit-grid" ]]; then
      combineTool.py \
      -M AsymptoticGrid \
      ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/vlq_asymptotic_grid.json \
      --redefineSignalPOIs r \
      --setParameterRanges r=0,1 \
      --setParameters r=1,lumi_scale=1 \
      -d ${datacarddir}/combined/cmb/ws.root \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name vlq_${sub_analysis}_full_cmb \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 \
      -v 1

  elif [[ $MODE == "run-impacts" ]]; then
      taskname="impacts_${TAG}_${sub_analysis}_mH${MASS}"

      combineTool.py -M Impacts -d ${datacarddir}/combined/cmb/ws.root \
      --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
      --doInitialFit --robustFit 1 \
      -m $MASS \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU  -v 0

      combineTool.py -M Impacts -d ${datacarddir}/combined/cmb/ws.root \
      --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 \
      --robustFit 1 --doFits \
      -m $MASS \
      --setParameters gU=0,lumi_scale=1 \
      --redefineSignalPOIs gU \
      --job-mode "SGE" \
      --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" \
      --task-name ${taskname} --merge 5

  elif [[ $MODE == "collect-impacts" ]]; then
      combineTool.py -M Impacts -d ${datacarddir}/combined/cmb/ws.root -m $MASS -o ${sub_analysis}_${MASS}_impacts.json --redefineSignalPOIs gU

  elif [[ $MODE == "ws-plot" ]]; then

      combineTool.py -M T2W -o "ws.root" \
      -P CombineHarvester.MSSMvsSMRun2Legacy.VLQ:VLQ \
      -i ${datacarddir}/201?/htt_*/ \
      --X-allow-no-signal \
      -m 1000 --parallel 4 | tee -a ${defaultdir}/logs/workspace_vlq.txt


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
      -n .combined-cmb.for_shape_unblinding \
      --there \
      -v 1

  elif [[ $MODE == "fit-for-scan" ]]; then

      combine -M MultiDimFit \
      -m $MASS \
      --setParameters gU=0.8885,lumi_scale=3.62  \
      -t -1 --freezeParameters gU \
      --redefineSignalPOIs gU \
      --X-rtd MINIMIZER_analytic \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
      --robustHesse 1 \
      -d ${datacarddir}/combined/cmb/ws.root \
      --points=30 --fastScan --algo=grid --rMin=0 --rMax=3 \
      -n "Scan"


  elif [[ $MODE == "lkld-plot" ]]; then
    plot1DScan.py higgsCombineScan.MultiDimFit.mH${MASS}.root --POI gU  --y-cut=8


  elif [[ $MODE == "prefit-plots" ]]; then
      freeze="MH=${MASS},gU=${gU}"
      for era in 2018; do
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

      for era in 2018; do
          bash plotting/plot_shapes_vlq.sh \
              ${era} \
              "${datacarddir}/combined/cmb/prefit_shapes_${freeze}.root" \
              "${datacarddir}/plots/prefit_shapes_$(echo ${freeze} | sed 's/=//g; s/\./p/g')/" \
              tt \
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

  elif [[ $MODE == "ccc-channel" ]]; then

      ulimit -s unlimited

      combine -M ChannelCompatibilityCheckRegexGroup \
      -d ${datacarddir}/combined/cmb/ws.root \
      -m $MASS \
      --setParameters gU=0,lumi_scale=1 \
      --setParameterRanges gU=-5,5 \
      --rMin=-5 \
      --redefineSignalPOIs gU \
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
      --setParameters gU=0,lumi_scale=1 \
      --setParameterRanges gU=-5,5 \
      --rMin=-5 \
      --redefineSignalPOIs gU \
      --X-rtd MINIMIZER_analytic \
      --X-rtd FITTER_NEW_CROSSING_ALGO \
      --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
      --robustFit 1 --stepSize 0.01 \
      -n .CCC.year \
      --saveFitResult \
      -g 2016 -g 2017 -g 2018 \
      -v 1


  elif [[ $MODE == "plot-ccc" ]]; then
      python plotting/plot_ccc.py \
      higgsCombine.CCC.combined-cmb.ChannelCompatibilityCheckRegexGroup.mH${MASS}.root \
      -o ChannelCompatibilityCheck_FitResults_mH${MASS} \
      -p gU \
      -r m5,5

  fi
done
