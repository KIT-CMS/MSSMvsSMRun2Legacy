# MSSMvsSMRun2Legacy
Analysis specific software for MSSM (with SM categories) to be used within CombineHarvester

# Setup software

```bash
wget https://raw.githubusercontent.com/KIT-CMS/MSSMvsSMRun2Legacy/master/scripts/checkout.sh
bash ./checkout.sh
```
After this, please put the datacard inputs into the [shapes](shapes) folder in the following structure:

`<year>/<channel>/htt_<category>.inputs-mssm-vs-sm-Run<year>-<variable>.root`

# Cut-based SM analysis

## Datacard creation

```bash
morph_parallel.py --output output --analysis sm --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt --variable mt_tot_puppi --parallel 5
morph_parallel.py --output output --analysis sm --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_categories.txt --variable m_sv_puppi --parallel 5

for era in 2016 2017 2018;
do
    mkdir -p output_sm/${era}/cmb/; rsync -av --progress output_sm/${era}/htt_*/*  output_sm/${era}/cmb/
done;
mkdir -p output_sm/combined/cmb/; rsync -av --progress output_sm/201?/htt_*/*  output_sm/combined/cmb/
```

## Workspace creation

**inclusive and stage 0:**

```bash
ulimit -s unlimited
combineTool.py -M T2W -i output_sm/{2016,2017,2018,combined}/cmb --parallel 4 -o inclusive.root
combineTool.py -M T2W -i output_sm/{2016,2017,2018,combined}/cmb --parallel 4 -o stage0.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"map=^.*/ggH125.?$:r_ggH[1,-9,11]"' --PO '"map=^.*/qqH125.?$:r_qqH[1,-9,11]"'
```

## Signal strength scans

**inclusive:**

```bash
combineTool.py -M MultiDimFit -d output_sm/combined/cmb/inclusive.root --there --algo grid --robustFit 1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --floatOtherPOIs 1 -t -1 --expectSignal 1 -v1 --alignEdges 1 --setParameterRanges r=0.85,1.15 --points 41 --split-points 1 -n .inclusive --job-mode condor --task-name r_1Dscan --dry-run

# after adapting the condor configuration, submit to batch system:
condor_submit condor_r_1Dscan.sub
```

**stage 0:**

```bash
# ggH
combineTool.py -M MultiDimFit -d output_sm/combined/cmb/stage0.root --there --algo grid --robustFit 1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --floatOtherPOIs 1 -t -1 --expectSignal 1 -v1 --setParameters r_ggH=1,r_qqH=1 --redefineSignalPOIs r_ggH --alignEdges 1 --setParameterRanges r_ggH=0.65,1.35:r_qqH=0.3,1.7 --points 41 --split-points 1 -n .stage0_r_ggH --job-mode condor --task-name r_ggH_1Dscan --dry-run

# qqH
combineTool.py -M MultiDimFit -d output_sm/combined/cmb/stage0.root --there --algo grid --robustFit 1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --floatOtherPOIs 1 -t -1 --expectSignal 1 -v1 --setParameters r_ggH=1,r_qqH=1 --redefineSignalPOIs r_qqH --alignEdges 1 --setParameterRanges r_qqH=0.65,1.35:r_ggH=0.3,1.7 --points 41 --split-points 1 -n .stage0_r_qqH --job-mode condor --task-name r_qqH_1Dscan --dry-run

# after adapting the condor configuration, submit to batch system:
condor_submit condor_r_ggH_1Dscan.sub
condor_submit condor_r_qqH_1Dscan.sub

# determining correlation ggH vs. qqH
combineTool.py -M FitDiagnostics -d output_sm/combined/cmb/stage0.root --there --robustHesse 1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -t -1 --expectSignal 1 -v1 --setParameters r_ggH=1,r_qqH=1 --setParameterRanges r_qqH=0.65,1.35:r_ggH=0.65,1.35 -n .stage0 --task-name stage0_correlation
```

**Collecting results & plotting likelihood scans:**

```bash
hadd output_sm/combined/cmb/higgsCombine.inclusive.MultiDimFit.mH120.root output_sm/combined/cmb/higgsCombine.inclusive.POINTS*.root
hadd output_sm/combined/cmb/higgsCombine.stage0_r_ggH.MultiDimFit.mH120.root output_sm/combined/cmb/higgsCombine.stage0_r_ggH.POINTS*.root
hadd output_sm/combined/cmb/higgsCombine.stage0_r_qqH.MultiDimFit.mH120.root output_sm/combined/cmb/higgsCombine.stage0_r_qqH.POINTS*.root

plot1DScan.py output_sm/combined/cmb/higgsCombine.inclusive.MultiDimFit.mH120.root -o r_1Dscan
plot1DScan.py output_sm/combined/cmb/higgsCombine.stage0_r_ggH.MultiDimFit.mH120.root -o r_ggH_1Dscan --POI r_ggH
plot1DScan.py output_sm/combined/cmb/higgsCombine.stage0_r_qqH.MultiDimFit.mH120.root -o r_qqH_1Dscan --POI r_qqH

print_correlation.py output_sm/combined/cmb/fitDiagnostics.stage0.root r_ggH r_qqH
```

## Prefit shapes

```bash
for card in output_MSSMvsSM_Run2_sm_nobtag_m_sv_puppi/*/cmb/combined.txt.cmb;
do
    PostFitShapesFromWorkspace -w ${card/combined.txt.cmb/stage0.root}  -o ${card/combined.txt.cmb/prefit_shapes.root} -d ${card}
done;
```

# Model-independent MSSM analysis

## Datacard creation

```bash
morph_parallel.py --output output --analysis mssm --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt --variable mt_tot_puppi --additional_arguments="--auto_rebin=1" --sm_gg_fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3.root --parallel 5
morph_parallel.py --output output --analysis mssm --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt --variable mt_tot_puppi --additional_arguments="--auto_rebin=1" --sm_gg_fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3.root --parallel 5

for era in 2016 2017 2018;
do
    mkdir -p output_mssm/${era}/cmb/; rsync -av --progress output_mssm/${era}/htt_*/*  output_mssm/${era}/cmb/
done;
mkdir -p output_mssm/combined/cmb/; rsync -av --progress output_mssm/201?/htt_*/*  output_mssm/combined/cmb/
```

## Workspace creation

```bash
ulimit -s unlimited
combineTool.py -M T2W -o "ws.root" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' -i output_mssm/{2016,2017,2018,combined}/cmb -m 110 --parallel 4
```

## Model-independent CLs 95% limits (asymptotic, SM Higgs in background hypothesis)

**bbH:**

```bash
combineTool.py -m "110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_bbH -d output_mssm/combined/cmb/ws.root --there -n ".bbH" --job-mode condor --dry-run --task-name bbH_full_combined --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1

# After adaption of condor_bbH_full_combined.sub, submit to batch system:
condor_submit condor_bbH_full_combined.sub
```

**ggH:**

```bash
combineTool.py -m "110,120,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_ggH -d output_mssm/combined/cmb/ws.root --there -n ".ggH" --job-mode condor --dry-run --task-name ggH_full_combined --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 -v 1

# After adaption of condor_ggH_full_combined.sub, submit to batch system:
condor_submit condor_ggH_full_combined.sub
```

**Collecting limits:**

```bash
for p in gg bb
do
    combineTool.py -M CollectLimits output_mssm/combined/cmb/higgsCombine.${p}H*.root --use-dirs -o output_mssm/combined/cmb/mssm_${p}H_combined.json
done
```

**Plotting limits:**

```bash
for p in gg bb
do
    plotMSSMLimits.py --cms-sub "Own Work" --title-right "137 fb^{-1} (13 TeV)" --process "${p}#phi" --y-axis-min 0.0001 --y-axis-max 1000.0 --show exp,obs output_mssm/combined/cmb/mssm_${p}H_combined_cmb.json  --output mssm_model-independent_combined_${p}H_cmb --logx --logy
done
```

## Plotting of prefit shapes

**Workspace creation:**
In order to create prefit shapes from the available datacards, separate workspaces for the different analysis categories have to be created. This can be done with the following command
```bash
ulimit -s unlimited
combineTool.py -M T2W -o "ws.root" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' -i output_mssm_classic/{2016,2017,2018}/htt_* -m 110 --parallel 5
```

**Prefit shapes:**
An exemplary command to extract prefit shapes from the created workspaces and datacards is
```bash
prefit_postfit_shapes_parallel.py --datacard_pattern "output_mssm_classic/201?/htt_*/combined.txt.cmb" --workspace_name ws.root --output_name prefit_shapes.root --freeze_arguments "--freeze MH=1200,r_ggH=0.1,r_bbH=0.1" --parallel 5

hadd output_mssm_classic/combined/cmb/prefit_shapes.root output_mssm_classic/{2016,2017,2018}/htt_*/prefit_shapes.root
```
where the freeze arguments can be chosen freely.

**Plotting of the created shapes:**
Plotting scripts are provided in the [`plotting`](https://github.com/KIT-CMS/MSSMvsSMRun2Legacy/tree/master/plotting) directory of this repository. The command to produce the plots of the prefit shapes created above is
```bash
plot_shapes_mssm_model_independent.sh $ERA output_mssm_classic/combined/cmb/prefit_shapes.root $OUTPUT_DIR et,mt,tt,em 1200
```
The last parameter given to the script is optional and will be the mass that is displayed in the legend of the plot. Currently the plotting is only possible for one experiment era a time. The output directory can be freely chosen and will be created if it is not yet present.

# Model-dependent MSSM analysis

In case of the calculation of model-dependent limits, several choices of signal models and categorization are available, denoted with an appropriate `--analysis`.
The full BSM (in this case, MSSM) signal model consists of three neutral Higgs bosons, the light scalar `h`, the heavy scalar `H`, and the heavy pseudoscalar `A`. The
production modes `ggphi` and `bbphi` are considered for all three Higgs bosons, whereas the `qqphi` production mode is accounted for only for the light scalar `h`.

The different analysis setups are listed in the following:

 * `mssm_vs_sm_classic`: using classic categorization (similar to HIG-17-020) with the **full** signal model for the BSM prediction. The corresponding hypothesis is then compared with the SM prediction.
 * `mssm_vs_sm_classic_h125`: same as above, but using for ggh the templates from the SM prediction, which are reweighted to the yield predicted by the MSSM scenario.
 * `mssm_vs_sm_heavy`: using classic categorization (similar to HIG-17-020), but using only the heavy Higgs boson predictions `H` and `A` added on top of the background and SM prediction. This hypothesis
is compared with the SM prediction.
 * `mssm_vs_sm`: using SM categorization in addition to classic categories, with the **full** signal model for the BSM prediction. The corresponding hypothesis is then compared with the SM prediction. The BSM signal modelling is dependent on the categories:
   * for SM categories, only ggh and qqh are taken into account
   * high mass no-btag categories contain only bbH, bbA and ggH, ggA
   * btag categories & common control regions contain the full BSM signal model
 * `mssm_vs_sm_h125`: same as above, but using for ggh the templates from the SM prediction, which are reweighted to the yield predicted by the MSSM scenario.
 * `mssm_vs_sm_CPV` : here the considered Higgs bosons are `H1` (SM-like), `H2` and `H3` instead of the usual `h` `H` and `A`. _IN PROGRESS_

## Datacard creation for `mssm_vs_sm_classic_h125`

```bash
morph_parallel.py --output output --analysis mssm_vs_sm_classic_h125 --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt --variable mt_tot_puppi --parallel 5
morph_parallel.py --output output --analysis mssm_vs_sm_classic_h125 --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt --variable mt_tot_puppi --parallel 5

mkdir -p output_mssm_vs_sm_classic_h125/combined/cmb/; rsync -av --progress output_mssm_vs_sm_classic_h125/201?/htt_*/*  output_mssm_vs_sm_classic_h125/combined/cmb/
```

## Datacard creation for `mssm_vs_sm_heavy`

```bash
morph_parallel.py --output output --analysis mssm_vs_sm_heavy --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt --variable mt_tot_puppi --parallel 5
morph_parallel.py --output output --analysis mssm_vs_sm_heavy --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt --variable mt_tot_puppi --parallel 5

mkdir -p output_mssm_vs_sm_heavy/combined/cmb/; rsync -av --progress output_mssm_vs_sm_heavy/201?/htt_*/*  output_mssm_vs_sm_heavy/combined/cmb/
```

## Datacard creation for `mssm_vs_sm_h125`

```bash
morph_parallel.py --output output --analysis mssm_vs_sm_h125 --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt --variable mt_tot_puppi --parallel 5
morph_parallel.py --output output --analysis mssm_vs_sm_h125 --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_new_categories.txt --variable mt_tot_puppi --parallel 5
morph_parallel.py --output output --analysis mssm_vs_sm_h125 --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/sm_categories.txt --variable m_sv_puppi --parallel 5

mkdir -p output_mssm_vs_sm_h125/combined/cmb/; rsync -av --progress output_mssm_vs_sm_h125/201?/htt_*/*  output_mssm_vs_sm_h125/combined/cmb/
```

## Datacard creation for `mssm_vs_sm_CPV`

```bash
morph_parallel.py --output output --analysis mssm_vs_sm_CPV --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/control_region_categories.txt --variable mt_tot_puppi --parallel 5 --sm_gg_fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3.root
morph_parallel.py --output output --analysis mssm_vs_sm_CPV --eras 2016,2017,2018 --category_list ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_classic_categories.txt --variable mt_tot_puppi --parallel 5 --sm_gg_fractions ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3.root

mkdir -p output_mssm_vs_sm_CPV/combined/cmb/; rsync -av --progress output_mssm_vs_sm_CPV/201?/htt_*/*  output_mssm_vs_sm_CPV/combined/cmb/
```

## Workspace creation (exemplary for `mssm_vs_sm_classic_h125`)

```bash
ulimit -s unlimited
combineTool.py -M T2W -o ws_mh125.root  -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM:MSSMvsSM --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ --PO modelFile=13,Run2017,mh125_13.root --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3_mssm_mode.root -i output_mssm_vs_sm_classic_h125/combined/cmb/ --PO minTemplateMass=110.0 --PO maxTemplateMass=3200.0
```
_NB_ : for the `mssm_vs_sm_CPV` analysis, use `mh1125_CPV_13.root` instead of `mh125_13.root` and `--PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3.root` instead of `--PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3_mssm_mode.root`. The model building file to use is `CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM_CPV` with the model named `MSSMvsSM_CPV`. You may also change `ws_mh125.root` to `ws_mH1125_CPV.root`. The command to run would then be
```
combineTool.py -M T2W -o ws_mH1125_CPV.root  -P CombineHarvester.MSSMvsSMRun2Legacy.MSSMvsSM_CPV:MSSMvsSM_CPV --PO filePrefix=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/ --PO modelFile=13,Run2017,mh1125_CPV_13.root --PO MSSM-NLO-Workspace=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_v3.root -i output_mssm_vs_sm_CPV/combined/cmb/ --PO minTemplateMass=110.0 --PO maxTemplateMass=3200.0
```

## Model-dependent CLs 95% limits for `mssm_vs_sm_classic_h125`

**Computing limits:**

```bash
mkdir -p calculation_mh125_mssm_vs_sm_classic_h125
cd calculation_mh125_mssm_vs_sm_classic_h125

ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_classic_h125/combined/cmb/ws_mh125.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_classic_h125_1' --dry-run --redefineSignalPOI x --setParameters r=1 --freezeParameters r -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01

# After adaption of each shell script and condor configuration matching mattern condor_mssm_mh125_mssm_vs_sm_classic_h125_1.{sh,sub}, submit to batch system:
condor_submit  condor_mssm_vs_sm_classic_h125_1.sub
```

**Collecting limits:**

Basically the same command as above, but with a different task name (in case some of the jobs have failed)

```bash
ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_classic_h125/combined/cmb/ws_mh125.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_classic_h125_2' --dry-run --redefineSignalPOI x --setParameters r=1 --freezeParameters r -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01
```

**Plotting limits:**

```bash
plotLimitGrid.py asymptotic_grid.root --scenario-label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)" --output mssm_mh125_mssm_vs_sm_classic_h125  --title-right="137 fb^{-1} (13 TeV)" --cms-sub="Own Work" --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/mh125_13.root --y-range 2.0,60.0 --x-title "m_{A} [GeV]"
```
## Model-dependent CLs 95% limits for `mssm_vs_sm_h125`

**Computing limits:**

```bash
mkdir -p calculation_mh125_mssm_vs_sm_h125
cd calculation_mh125_mssm_vs_sm_h125

ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_h125/combined/cmb/ws_mh125.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_h125_1' --dry-run --redefineSignalPOI x --setParameters r=1 --freezeParameters r -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01

# After adaption of each shell script and condor configuration matching mattern condor_mssm_mh125_mssm_vs_sm_h125_1.{sh,sub}, submit to batch system:
condor_submit  condor_mssm_vs_sm_h125_1.sub
```

**Collecting limits:**

Basically the same command as above, but with a different task name (in case some of the jobs have failed)

```bash
ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_h125/combined/cmb/ws_mh125.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_h125_2' --dry-run --redefineSignalPOI x --setParameters r=1 --freezeParameters r -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01
```

**Plotting limits:**

```bash
plotLimitGrid.py asymptotic_grid.root --scenario-label="M_{h}^{125} scenario (h,H,A#rightarrow#tau#tau)" --output mssm_mh125_mssm_vs_sm_h125  --title-right="137 fb^{-1} (13 TeV)" --cms-sub="Own Work" --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/mh125_13_fixed.root --y-range 2.0,60.0 --x-title "m_{A} [GeV]"
```

## Model-dependent CLs 95% limits for `mssm_vs_sm_heavy`

**Computing limits:**

```bash
mkdir -p calculation_mh125_mssm_vs_sm_heavy
cd calculation_mh125_mssm_vs_sm_heavy

ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_heavy/combined/cmb/ws_mh125.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_heavy_1' --dry-run -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01

# After adaption of each shell script and condor configuration matching mattern condor_mssm_mh125_mssm_vs_sm_heavy_1.{sh,sub}, submit to batch system:
condor_submit  condor_mssm_vs_sm_heavy_1.sub
```

**Collecting limits:**

Basically the same command as above, but with a different task name (in case some of the jobs have failed)

```bash
ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_mh125.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_heavy/combined/cmb/ws_mh125.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_heavy_2' --dry-run --freezeParameters r -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01
```

**Plotting limits:**

```bash
plotLimitGrid.py asymptotic_grid.root --scenario-label="M_{h}^{125} scenario (H,A#rightarrow#tau#tau)" --output mssm_mh125_mssm_vs_sm_heavy  --title-right="137 fb^{-1} (13 TeV)" --cms-sub="Own Work" --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/mh125_13_fixed.root --y-range 2.0,60.0 --x-title "m_{A} [GeV]"
```

## Model-dependent CLs 95% limits for `mssm_vs_sm_CPV`

**Computing limits:**

```bash
mkdir -p calculation_mh125_mssm_vs_sm_CPV
cd calculation_mh125_mssm_vs_sm_CPV

ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_CPV.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_CPV/combined/cmb/ws_mH1125_CPV.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_CPV_1' --dry-run -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01 --redefineSignalPOI x --setParameters r=1 --freezeParameters r

# After adaption of each shell script and condor configuration matching mattern condor_mssm_mh125_mssm_vs_sm_CPV_1.{sh,sub}, submit to batch system:
condor_submit  condor_mssm_mh125_mssm_vs_sm_CPV_1.sub
```

**Collecting limits:**

Basically the same command as above, but with a different task name (in case some of the jobs have failed)

```bash
ulimit -s unlimited
combineTool.py -M AsymptoticGrid ${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_asymptotic_grid_CPV.json -d ${CMSSW_BASE}/src/output_mssm_vs_sm_CPV/combined/cmb/ws_mH1125_CPV.root --job-mode 'condor' --task-name 'mssm_mh125_mssm_vs_sm_CPV_2' --dry-run -v1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerTolerance 0.01 --redefineSignalPOI x --setParameters r=1 --freezeParameters r
```

**Plotting limits:**

```bash
plotLimitGrid.py asymptotic_grid.root --scenario-label="M_{H_{1}}^{125} CPV scenario (H_{1},H_{2},H_{3}#rightarrow#tau#tau)" --output mssm_mh125_mssm_vs_sm_CPV  --title-right="137 fb^{-1} (13 TeV)" --cms-sub="Own Work" --contours="exp-2,exp-1,exp0,exp+1,exp+2,obs" --model_file=${CMSSW_BASE}/src/CombineHarvester/MSSMvsSMRun2Legacy/data/mh1125_CPV_13.root --y-range 1.0,20.0 --x-title "m_{H^{#pm}} [GeV]" --mass_histogram m_H1
```

## Adaptions needed for job submission via htcondor

The following line in file [CombineToolBase.py](https://github.com/KIT-CMS/CombineHarvester/blob/master/CombineTools/python/combine/CombineToolBase.py#L33) from the CombineHarvester package needs to be
updated to be able to run on htcondor batch systems with the following configuration snippet:

**ETP batch system**

```python
Requirements = ( (Target.ProvidesCPU == True) && (TARGET.ProvidesEKPResources == True ) )
universe = docker
docker_image = mschnepf/slc7-condocker
+RequestWalltime = 10800
+ExperimentalJob = True
RequestMemory = 2000
RequestCpus = 1
accounting_group = cms.higgs
```

**No additions needed for CERN or NAF batch systems**
