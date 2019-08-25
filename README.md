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
