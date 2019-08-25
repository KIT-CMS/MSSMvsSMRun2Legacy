# MSSMvsSMRun2Legacy
Analysis specific software for MSSM (with SM categories) to be used within CombineHarvester

# Cut-based SM analysis

## Datacard creation

```bash
MorphingMSSMvsSM --era=2017 --auto_rebin=1 --binomial_bbb=1 --variable=m_sv_puppi --categories=sm_nobtag
```

## Workspace creation

**inclusive and stage 0:**

```bash
combineTool.py -M T2W -i output_MSSMvsSM_Run2_sm_nobtag_m_sv_puppi/*/* --parallel 10 -o inclusive.root
combineTool.py -M T2W -i output_MSSMvsSM_Run2_sm_nobtag_m_sv_puppi/*/* --parallel 10 -o stage0.root -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"map=^.*/ggH125.?$:r_ggH[1,-9,11]"' --PO '"map=^.*/qqH125.?$:r_qqH[1,-9,11]"'
```

## Signal strength fits

**inclusive:**

```bash
combineTool.py -M MultiDimFit -d output_MSSMvsSM_Run2_sm_nobtag_m_sv_puppi/*/*/inclusive.root --algo singles --robustFit 1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --floatOtherPOIs 1 -t -1 --expectSignal 1 -n .inclusive -v1 --there --parallel 10 --setParameterRanges r=-3.0,5.0
```

**stage 0:**

```bash
combineTool.py -M MultiDimFit -d output_MSSMvsSM_Run2_sm_nobtag_m_sv_puppi/*/*/stage0.root --algo singles --robustFit 1 --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --floatOtherPOIs 1 -t -1 --expectSignal 1 -n .stage0 -v1 --there --parallel 10
```

**Printing out fit results:**

```bash
for result in output_MSSMvsSM_Run2_sm_nobtag_m_sv_puppi/*/*/higgsCombine.*.MultiDimFit*.root;
do
    print_fitresult.py ${result};
    echo;
done;
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
MorphingMSSMvsSM --era=2017 --auto_rebin=1 --binomial_bbb=1 --variable=mt_tot_puppi --categories=mssm
```

## Workspace creation

```bash
combineTool.py -M T2W -o "ws.root" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO '"map=^.*/ggh_(i|t|b).?$:r_ggH[0,0,200]"' --PO '"map=^.*/bbh$:r_bbH[0,0,200]"' -i output_MSSMvsSM_Run2_mssm_mt_tot_puppi/Run2017/cmb/
```

## Model-independent CL 95% limits (asymptotic)

**bbH:**

```bash
combineTool.py -m "100,110,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_bbH -d output_MSSMvsSM_Run2_mssm_mt_tot_puppi/Run2017/cmb/ws.root --there -n ".bbH" --job-mode condor --dry-run --task-name bbH_full

# After adaption of condor_bbH_full.sh and condor_bbH_full.sub; submit to batch system:
condor_submit condor_bbH_full.sub
```

**ggH**

```bash
combineTool.py -m "100,110,120,125,130,140,160,180,200,250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200" -M AsymptoticLimits --rAbsAcc 0 --rRelAcc 0.0005 --boundlist CombineHarvester/MSSMvsSMRun2Legacy/input/mssm_boundaries.json --setParameters r_ggH=0,r_bbH=0 --redefineSignalPOIs r_ggH -d output_MSSMvsSM_Run2_mssm_mt_tot_puppi/Run2017/cmb/ws.root --there -n ".ggH" --job-mode condor --dry-run --task-name ggH_full

# After adaption of condor_ggH_full.sh and condor_ggH_full.sub; submit to batch system:
condor_submit condor_ggH_full.sub
```

**Collecting limits:**

```bash
for p in ggH bbH
do
    combineTool.py -M CollectLimits output_MSSMvsSM_Run2_mssm_mt_tot_puppi/Run2017/cmb/higgsCombine.${p}*.root --use-dirs -o output_MSSMvsSM_Run2_mssm_m_sv_puppi/Run2017/cmb/mssm_${p}.json
done
```

**Plotting limits:**

```bash
for p in gg bb
do
    plotMSSMLimits.py --cms-sub "Private Work" --title-right "41.5 fb^{-1} (2017, 13 TeV)" --process "${p}#phi" --y-axis-min 0.001 --y-axis-max 1000.0 --show exp,obs output_MSSMvsSM_Run2_mssm_mt_tot_puppi/Run2017/cmb/mssm_${p}H_cmb.json  --output mssm_m_sv_${p}H_cmb --logx --logy
done
```
