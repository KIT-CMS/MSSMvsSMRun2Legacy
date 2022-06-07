#!/bin/bash

ulimit -s unlimited

TAG=$1
MODE=$2
if [[ $TAG == "auto" ]]; then
    TAG="cmb_ind_lowmass"
fi

xrange () {
    [ "$#" != 2 ] && exit
    case $1 in
        "60")
            echo 8.
            ;;
        "80")
            if [[ $2 == 1 ]]; then echo 5.; else echo 3.; fi
            ;;
        "95")
            if [[ $2 == 1 ]]; then echo 25.; else echo 30.; fi
            ;;
        "100")
            if [[ $2 == 1 ]]; then echo 15.; else echo 20.; fi
            ;;
        "120")
            if [[ $2 == 1 ]]; then echo 3.; else echo 6.; fi
            ;;
        "125")
            if [[ $2 == 1 ]]; then echo 2.; else echo 4.5; fi
            ;;
        "130")
            if [[ $2 == 1 ]]; then echo 1.2; else echo 4.; fi
            ;;
        "140")
            if [[ $2 == 1 ]]; then echo 0.8; else echo 2.5; fi
            ;;
        "160")
            if [[ $2 == 1 ]]; then echo 0.5; else echo 1.; fi
            ;;
        "180")
            if [[ $2 == 1 ]]; then echo 0.4; else echo 0.3; fi
            ;;
        "200")
            if [[ $2 == 1 ]]; then echo 0.3; else echo 0.15; fi
            ;;
        *)
            exit 1
            ;;
    esac
}

yrange () {
    [ "$#" != 2 ] && exit
    case $1 in
        "60")
            echo 50.
            ;;
        "80")
            if [[ $2 == 1 ]]; then echo 22.; else echo 16.; fi
            ;;
        "95")
            if [[ $2 == 1 ]]; then echo 6.; else echo 5.; fi
            ;;
        "100")
            echo 4.
            ;;
        "120")
            if [[ $2 == 1 ]]; then echo 1.; else echo 1.5; fi
            ;;
        "125")
            if [[ $2 == 1 ]]; then echo 0.8; else echo 1.2; fi
            ;;
        "130")
            if [[ $2 == 1 ]]; then echo 0.6; else echo 1.0; fi
            ;;
        "140")
            if [[ $2 == 1 ]]; then echo 0.5; else echo 0.8; fi
            ;;
        "160")
            echo 0.4
            ;;
        "180")
            if [[ $2 == 1 ]]; then echo 0.3; else echo 0.25; fi
            ;;
        "200")
            if [[ $2 == 1 ]]; then echo 0.35; else echo 0.25; fi
            ;;
        *)
            exit 1
            ;;
    esac
}

defaultdir=analysis_2022_02_28/$TAG
analysis="bsm-model-indep"
hSM_treatment="hSM-in-bg"
# hSM_treatment="no-hSM-in-bg"
# categorization="lowmass"
# sm_like_hists="sm125"
real_data=1
[[ ! -d ${defaultdir} ]] && mkdir -p ${defaultdir}
[[ ! -d ${defaultdir}/logs ]] && mkdir -p ${defaultdir}/logs
[[ ! -d ${defaultdir}/limits/condor ]] && mkdir -p ${defaultdir}/limits/condor
[[ ! -d ${defaultdir}/limits_ind/condor ]] && mkdir -p ${defaultdir}/limits_ind/condor
defaultdir=$(readlink -f analysis_2022_02_28/$TAG)
datacarddir=${defaultdir}/datacards_${analysis}
identifier_toy_submit=$(date +%Y_%m_%d)
[[ -z $3 ]] || identifier_toy_submit=$3

case "$MODE" in

    "initial")
        morph_parallel.py \
            --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis "none" \
            --hSM-treatment ${hSM_treatment} \
            --categorization="lowmass" \
            --sm-like-hists="sm125" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_2d_to_1d.txt \
            --variable "m_sv_VS_pt_tt_splitpT" \
            --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=${real_data}"
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
            --parallel 5 |& tee ${defaultdir}/logs/morph_mssm_log.txt

        # take btag cats from m_sv binned one
        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis "none" \
            --hSM-treatment ${hSM_treatment} \
            --categorization="lowmass" \
            --sm-like-hists="sm125" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_1d_btag.txt \
            --variable "m_sv_puppi" \
            --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=${real_data}" \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
            --parallel 5 |& tee -a ${defaultdir}/logs/morph_mssm_log.txt


        morph_parallel.py --output ${defaultdir}/datacards \
            --analysis ${analysis} \
            --sub-analysis "none" \
            --hSM-treatment ${hSM_treatment} \
            --categorization="lowmass" \
            --sm-like-hists="sm125" \
            --eras 2016,2017,2018 \
            --category-list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories_cr.txt \
            --variable "mt_tot_puppi" \
            --additional-arguments="--auto_rebin=1 --manual_rebin=1 --real_data=${real_data}" \
            --sm-gg-fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root \
            --parallel 5 |& tee -a ${defaultdir}/logs/morph_mssm_log.txt

        #################
        # combining output
        #################
        mkdir -p ${datacarddir}/combined/cmb/
        rsync -av --progress ${datacarddir}/201?/htt_*/* ${datacarddir}/combined/cmb/ |& tee ${defaultdir}/logs/copy_datacards.txt
        for year in 2016 2017 2018; do
            mkdir -p ${datacarddir}/${year}/cmb/
            rsync -av --progress ${datacarddir}/${year}/htt_*/* ${datacarddir}/${year}/cmb/ |& tee -a ${defaultdir}/logs/copy_datacards.txt
        done

        # Check if all datacards that should be present have been written
        # Number of categories: et+mt+tt+em
        EXPECTED=$(((5+5+5+11)*3))
        if [[ $(ls ${datacarddir}/combined/cmb/*.txt | wc -l) != $EXPECTED ]]; then
            echo -e "\033[0;31m[ERROR]\033[0m Not all datacards have been created or written. Please check the logs..."
            echo "Expected ${EXPECTED} datacards written but found only $(ls ${datacarddir}/combined/cmb/ | wc -l) in the combined directory."
        fi
        ;;

    "ws")
        #################
        # workspace creation
        #################
        combineTool.py -M T2W -o "ws.root" \
            -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
            --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' \
            --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' \
            --PO '"map=^.*/qqX$:r_qqX[0]"' \
            --PO '"map=^.*/ggX_(i|t|b).?$:r_ggX[0,0,200]"' \
            -i ${datacarddir}/combined/cmb \
            -m 95 --parallel 8
      ;;

    "prepare-ggH-bbH-scan")
        [[ ! -d ${defaultdir}/ggH_bbH_scan_ind/condor ]] && mkdir -p ${defaultdir}/ggH_bbH_scan_ind/condor
        cd ${defaultdir}/ggH_bbH_scan_ind/condor
        # Run 2D likelihood scans for r_ggH and r_bbH
        combineTool.py -M MultiDimFit \
            --algo grid --points 40000 --alignEdges 1 --split-points 50 \
            --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_ggH_bbH_2D_boundaries.json \
            -m "60,80,95,100,120,125,130,140,160,180,200" \
            --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH,r_bbH --freezeParameters r_qqX,r_ggX \
            -d ${datacarddir}/combined/cmb/ws.root \
            --job-mode condor --dry-run --task-name ggH_bbH_likelihood_scan \
            --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
            --cminFallbackAlgo Minuit2,Migrad,0:0.1 \
            -n ".ggH-bbH" \
            -v 1
            # --algo grid --points 100 --alignEdges 1 --split-points 10 \
            # -m "60,100,125,160" \
        ;;

    "submit-ggH-bbH-scan")
        cd ${defaultdir}/ggH_bbH_scan_ind/condor
        condor_submit condor_ggH_bbH_likelihood_scan.sub
        ;;

    "submit-ggH-bbH-scan-gc")
        gcworkdir=${defaultdir}/ggH_bbH_scan_ind/gc_condor_${identifier_toy_submit}
        mkdir -p ${gcworkdir}
        python scripts/build_gc_job.py \
            --combine-script ${defaultdir}/ggH_bbH_scan_ind/condor/condor_ggH_bbH_likelihood_scan.sh \
            --workspace ${datacarddir}/combined/cmb/ws.root \
            --workdir ${gcworkdir} \
            --tag ggH_bbH_likelihood_scan \
            --se-path /storage/gridka-nrg/$(whoami)/gc_storage/combine/ggH_bbH_likelihood_scan_lowmass_${identifier_toy_submit}

        echo "Submit with ${CMSSW_BASE}/src/grid-control/go.py ${gcworkdir}/ggH_bbH_likelihood_scan.conf -Gc -m 3"
        # ${CMSSW_BASE}/src/grid-control/go.py ${gcworkdir}/ggH_bbH_likelihood_scan_${identifier_toy_submit}.conf -G -m 3
        ;;

    "copy-results-ggH-bbH-scan")
        rsync -avhP /storage/gridka-nrg/$(whoami)/gc_storage/combine/ggH_bbH_likelihood_scan_lowmass_${identifier_toy_submit}/output/ ${defaultdir}/ggH_bbH_scan_ind/condor
        ;;

    "collect-ggH-bbH-scan")
        cd ${defaultdir}/ggH_bbH_scan_ind/
        smexp=""
        # for mass in 60 100 125 160; do
        for mass in 60 80 95 100 120 125 130 140 160 180 200; do
            [[ "$TAG" =~ NohSMinBackground$ ]] && smexp="--sm-exp ${datacarddir}/combined/cmb/higgsCombine.2D.SM1.bestfit.MultiDimFit.mH${mass}.root"
            for cmssub in "" "Supplementary"; do
                cmssubadd="_CMS"
                [[ "$cmssub" == "Preliminary" ]] && cmssubadd=""
                [[ "$cmssub" == "Supplementary" ]] && cmssubadd="_Supplementary"
                isBkg=1
                [[ "$TAG" =~ NohSMinBackground ]] && isBkg=0
                for int in 0 1; do
                    if [[ $int -eq 1 ]]; then
                        int_arg="--interpolate-missing"
                        int_sub="_interpolated"
                    else
                        int_arg=""
                        int_sub=""
                    fi
                    python ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/plotting/plotMultiDimFit.py \
                        --title-right="138 fb^{-1} (13 TeV)" \
                        --cms-sub=${cmssub} \
                        --mass $mass \
                        -o 2D_limit_mH${mass}${cmssubadd}${int_sub} \
                        ${smexp} \
                        --debug-output test_mH${mass}${int_sub}.root \
                        ${int_arg} \
                        --x-title "#sigma#font[42]{(gg#phi)}#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)" \
                        --y-title "#sigma#font[42]{(bb#phi)}#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)" \
                        --x-axis-max $(xrange $mass $isBkg) --y-axis-max $(yrange $mass $isBkg) \
                        condor/higgsCombine.ggH-bbH.POINTS.*.mH${mass}.root
                        # --likelihood-database \
                        # --add-3sigma-contour \
                done
            done
        done
        ;;

    "ggH-bbH-scan-SM-expectation")
        # TODO: This option only makes sense for the with hSM not in the background
        # Create asimov dataset for SM-expectation.
        combineTool.py -M MultiDimFit \
            --algo none \
            -m "125" \
            --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH,r_bbH --freezeParameters r_qqX,r_ggX \
            --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_ggH_bbH_2D_boundaries.json \
            -d $(dirname $defaultdir)/cmb_ind_lowmass_hSMinBackground/datacards_bsm-model-indep/combined/cmb/ws.root \
            --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
            -t -1 --saveToys \
            --there -n ".2D.ToyDataset.SM1" \
            -v 1

        # Copy it to correct location
        rsync -av --progress $(dirname $defaultdir)/cmb_ind_lowmass_hSMinBackground/datacards_bsm-model-indep/combined/cmb/higgsCombine.2D.ToyDataset.SM1.MultiDimFit.mH125.123456.root ${datacarddir}/combined/cmb/higgsCombine.2D.ToyDataset.SM1.MultiDimFit.mH125.123456.root

        # Run fits on this asimov dataset
        pushd ${defaultdir}/ggH_bbH_scan_ind/
        combineTool.py -M MultiDimFit \
            --algo none \
            -m "60,80,95,100,120,125,130,140,160,180,200" \
            --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_ggH_bbH_2D_boundaries_NohSMinBackground.json \
            --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 \
            --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH,r_bbH --freezeParameters r_qqX,r_ggX \
            -d ${datacarddir}/combined/cmb/ws.root \
            -t -1 --toysFile higgsCombine.2D.ToyDataset.SM1.MultiDimFit.mH125.123456.root \
            --there -n ".2D.SM1.bestfit"
            # --job-mode condor --dry-run --task-name ggH_bbH_likelihood_SM1 --merge 3 \
        ;;

    *)
        echo "[ERROR] Given mode $MODE not known..."
        exit 42
        ;;
esac
