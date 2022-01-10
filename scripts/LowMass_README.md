## Run limits and significances

```bash

# run morphing, T2W, and submit jobs for limits and significances

python run_model_independent_limits_2D_splitpTbins.py -o Dec20_2D

# collect limits into a json

for p in gg bb; do combineTool.py -M CollectLimits model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/higgsCombine.${p}H.v2.AsymptoticLimits.mH*.root --use-dirs -o model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/mssm_${p}H_combined.json; done

# plot limits

for p in gg bb; do plotMSSMLimits.py --cms-sub "" --title-right "138 fb^{-1} (13 TeV)" --process "${p}#phi" --y-axis-min 0.01 --y-axis-max 1000.0 --show exp,obs model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/mssm_${p}H_combined_cmb.json  --output mssm_model-independent_combined_${p}H_cmb  --logy; done

# plot significances
python scripts/significance_plot.py
```

## Run best fits for postfit plots

**With signal = 100 GeV:**

```bash
combineTool.py -m 100 -M MultiDimFit --saveFitResult  --freezeParameters r_qqX,r_ggX,r_bbH --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root --there -n "ggH.m100.bestfit.robustHesse"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --robustHesse=1
```

**Background only:**

```bash
combineTool.py -m 100 -M MultiDimFit --saveFitResult  --freezeParameters r_qqX,r_ggX,r_bbH,r_ggH --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root --there -n "ggH.bkgOnly.bestfit.robustHesse"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --robustHesse=1
```

## make post-fit / pre-fit plots

```bash

# run T2W on restore binning directory

combineTool.py -M T2W -i model_independent_limits/Dec20_2D_all_all_bsm-model-indep/restore_binning/

# make postfit plots with S+B fit

PostFitShapesFromWorkspace -w model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/restore_binning/combined.txt.cmb --fitresult model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/multidimfitggH.m100.bestfit.robustHesse.root:fit_mdf -o shapes_postfit_m100.root --skip-prefit=true   --mass 100 --postfit --freeze MH=100

# make postfit plots with B-only

PostFitShapesFromWorkspace -w model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/restore_binning/combined.txt.cmb --fitresult model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/multidimfitggH.bkgOnly.bestfit.robustHesse.root:fit_mdf -o shapes_postfit_bkgOnly.root --skip-prefit=true --mass 100 --postfit --freeze r_bbH=1,r_ggH=1,MH=100

# make prefit plots

PostFitShapesFromWorkspace -w model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/restore_binning/combined.txt.cmb  -o shapes_prefit.root --mass 100 --freeze r_bbH=1,r_ggH=1,MH=100

```

## Run approx impacts

**For ggH:**

```bash
# run fits

combineTool.py -M Impacts  --doFits --robustFit=1 --approx robust -m 100  --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root  -n ".ggH.impacts.dec20"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01

# collect results into a json

combineTool.py -M Impacts --approx robust -m 100  --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root  -n ".ggH.impacts.dec18" -o impacts_ggH_dec20.json

# make plot

plotImpacts.py -i impacts_ggH_dec20.json -o impacts_ggH_dec20 
```

**For bbH:**

```bash
# run fits

combineTool.py -M Impacts  --doFits --robustFit=1 --approx robust -m 100  --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_bbH -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root  -n ".bbH.impacts.dec20"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01

# collect results into a json

combineTool.py -M Impacts --approx robust -m 100  --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_bbH -d model_independent_limits/Dec20_2D_all_all_bsm-model-indep/combined/cmb/ws.root  -n ".bbH.impacts.dec18" -o impacts_bbH_dec20.json

# make plot

plotImpacts.py -i impacts_bbH_dec20.json -o impacts_bbH_dec20 
```

# full impacts

```bash

# run initial fit

combineTool.py -M Impacts  --doInitialFit --robustFit=1 -m 100  --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d model_independent_limits/Jan03_newresuncerts_all_all_bsm-model-indep/combined/cmb/ws.root  -n ".ggH.impacts_full.jan03"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01

# run seperate fits

combineTool.py -M Impacts  --doFits --robustFit=1 -m 100  --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d model_independent_limits/Jan03_newresuncerts_all_all_bsm-model-indep/combined/cmb/ws.root  -n ".ggH.impacts_full.jan03"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --named CMS_res_t > impacts_full_fits.out

# to run all parameters remove --named option

# collect results into json

combineTool.py -M Impacts   --robustFit=1 -m 100  --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --redefineSignalPOIs r_ggH -d model_independent_limits/Jan03_newresuncerts_all_all_bsm-model-indep/combined/cmb/ws.root  -n ".ggH.impacts_full.jan03"  --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -o impacts_full.json

# make plot

plotImpacts.py -i impacts_full.json -o impacts_full

```

## Channel compatability checks:

**By channel:**

```bash
MASS=100
datacarddir="model_independent_limits/Dec20_2D_all_all_bsm-model-indep/"

combine -M ChannelCompatibilityCheck \
-d ${datacarddir}/combined/cmb/ws.root \
-m $MASS \
--setParameters r_ggH=0,r_bbH=0,r_ggX=0,r_qqX=0  \
--freezeParameters r_ggX,r_qqX  \
--setParameterRange r_ggH=-20,20:r_bbH=0,100:CMS_htt_ttbarShape=-1.,1. \
--rMin=-20 \
--redefineSignalPOIs r_ggH \
--X-rtd MINIMIZER_analytic \
--X-rtd FITTER_NEW_CROSSING_ALGO \
--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
--robustFit 1 --stepSize 0.01 \
-n .CCC.channel \
--saveFitResult \
-g htt_em -g htt_et -g htt_mt -g htt_tt \
-v 2
```

**By era:**

```bash
MASS=100
datacarddir="model_independent_limits/Dec20_2D_all_all_bsm-model-indep/"

#run fits

combine -M ChannelCompatibilityCheck \
-d ${datacarddir}/combined/cmb/ws.root \
-m $MASS \
--setParameters r_ggH=0,r_bbH=0,r_ggX=0,r_qqX=0  \
--freezeParameters r_ggX,r_qqX  \
--setParameterRange r_ggH=-20,20:r_bbH=0,100:CMS_htt_ttbarShape=-1.,1. \
--rMin=-20 \
--redefineSignalPOIs r_ggH \
--X-rtd MINIMIZER_analytic \
--X-rtd FITTER_NEW_CROSSING_ALGO \
--cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.1 \
--robustFit 1 --stepSize 0.01 \
-n .CCC.year \
--saveFitResult \
-g 2016 -g 2017 -g 2018 \
-v 2

# make plot

python plotting/plot_ccc.py higgsCombine.CCC.year.ChannelCompatibilityCheck.mH${MASS}.root -p r_ggH -o ChannelCompatibilityCheck_FitResults_mH${MASS}.channel -r m12,15
```

## run GOF tests

```bash
# run fits on data and on toys for 2 scenarios:
# b-only (corresponds to M=95 below)
# and S+B fit with M=100 GeV (M=100 below)

DIR="Jan04_newemQCD_all_all_bsm-model-indep"

algo="saturated"

combineTool.py -m "95" -M GoodnessOfFit --algorithm ${algo} --boundlist input/mssm_boundaries.json  --there -d model_independent_limits/${DIR}/combined/cmb/ws.root -n ".${algo}" --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --freezeParameters r_qqX,r_ggX,r_bbH,r_ggH  --toysFreq --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" --task-name gof_${algo}_cmb_obs_${DIR}

combineTool.py -m "95" -M GoodnessOfFit --algorithm ${algo} --boundlist input/mssm_boundaries.json  --there -d model_independent_limits/${DIR}/combined/cmb/ws.root -n ".${algo}.toys" --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --freezeParameters r_qqX,r_ggX,r_bbH,r_ggH  --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 1 -s 0:499:1 --task-name gof_${algo}_cmb_toys_${DIR}

combineTool.py -m "100" -M GoodnessOfFit --algorithm ${algo} --boundlist input/mssm_boundaries.json  --there -d model_independent_limits/${DIR}/combined/cmb/ws.root -n ".${algo}" --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --freezeParameters r_qqX,r_ggX  --toysFreq --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" --task-name gof_${algo}_cmb_obs_${DIR}_sig

combineTool.py -m "100" -M GoodnessOfFit --algorithm ${algo} --boundlist input/mssm_boundaries.json  --there -d model_independent_limits/${DIR}/combined/cmb/ws.root -n ".${algo}.toys" --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0 --freezeParameters r_qqX,r_ggX  --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" -t 1 -s 0:499:1 --task-name gof_${algo}_cmb_toys_${DIR}_sig

# collect results

for m in 95 100; do combineTool.py -M CollectGoodnessOfFit --input model_independent_limits/${DIR}/combined/cmb/higgsCombine.${algo}.GoodnessOfFit.mH${m}.root model_independent_limits/${DIR}/combined/cmb/higgsCombine.${algo}.toys.GoodnessOfFit.mH${m}.*.root --there -o model_independent_limits/${DIR}/combined/cmb/gof_${algo}_m${m}.json; done

# make plots

```
## run scan of a parameter

```bash
# run scans
combineTool.py -m "100" -M MultiDimFit --boundlist input/mssm_boundaries.json --freezeParameters r_qqX,r_ggX --setParameters r_ggH=0,r_bbH=0,r_qqX=0,r_ggX=0,CMS_res_t=0 --redefineSignalPOIs CMS_res_t --setParameterRanges CMS_res_t=-2,2  -d model_independent_limits/Dec31_newsysts_tresonly_all_all_bsm-model-indep/combined/cmb/ws.root --there -n ".ggH.tres" --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01  --job-mode 'SGE' --prefix-file ic --sub-opts "-q hep.q -l h_rt=3:0:0" --split-points 1 --points 41 --algo grid --saveNLL --alignEdges=1

# hadd outputs

hadd -f model_independent_limits/Jan03_newresuncerts_all_all_bsm-model-indep/combined/cmb/higgsCombine.ggH.tres.MultiDimFit.mH100.root model_independent_limits/Jan03_newresuncerts_all_all_bsm-model-indep/combined/cmb/higgsCombine.ggH.tres.POINTS.*.root

# make plots

python scripts/plot1DScan.py --obs model_independent_limits/Jan03_newresuncerts_all_all_bsm-model-indep/combined/cmb/higgsCombine.ggH.tres.MultiDimFit.mH100.root --POI CMS_res_t --exp model_independent_limits/Jan03_newresuncerts_all_all_bsm-model-indep/combined/cmb/higgsCombine.ggH.tres.MultiDimFit.mH100.root --y-max=20 --y-cut=20 --output scan_tres 

```
