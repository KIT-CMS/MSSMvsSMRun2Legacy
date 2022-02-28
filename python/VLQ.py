from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import numpy as np
import ROOT

class VLQ(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        self.interference = True
        self.grid = False
        self.mu = False

    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("no_interference"):
                self.interference = False
            if po.startswith("grid"):
                self.grid = True
            if po.startswith("mu"):
                self.mu = True

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        poiNames = []

        if self.interference:
          scale = 1
        else:
          scale = 0

        if self.grid:
          self.modelBuilder.doVar("r[1,0,20]")
          poiNames.append('r')

        self.modelBuilder.doVar('gU[0,0,10]')
        poiNames.append('gU')

        #self.modelBuilder.doVar('r_ggH[0,0,10]')
        #poiNames.append('r_ggH')


        if self.mu:
          self.modelBuilder.doVar('mu[0,0,10]')
          poiNames.append('mu')

        self.modelBuilder.doSet('POI', ','.join(poiNames))

        #reweight_file = "/vols/cms/gu18/CH_unblinding/CMSSW_10_2_25/src/CombineHarvester/MSSMvsSMRun2Legacy/data/higgs_pt_reweighting_fullRun2_v2.root"
        #rf = ROOT.TFile(reweight_file)
        #w = rf.Get("w")
        #w.var("Yb_MSSM_H").setVal(1)
        #w.var("Yt_MSSM_H").setVal(1)
        #w.var("mH").setVal(1200)
        #Yt_scale = w.function("ggH_t_MSSM_frac").getVal()
        #Yb_scale = w.function("ggH_b_MSSM_frac").getVal()
        #Yi_scale = w.function("ggH_i_MSSM_frac").getVal()
        #self.modelBuilder.factory_('expr::ggH_t_scale("{}*@0", r_ggH)'.format(Yt_scale))
        #self.modelBuilder.factory_('expr::ggH_b_scale("{}*@0", r_ggH)'.format(Yb_scale))
        #self.modelBuilder.factory_('expr::ggH_i_scale("{}*@0", r_ggH)'.format(Yi_scale))

        if not self.mu:
          if self.grid:
            self.modelBuilder.factory_('expr::gU2("-{}*@0*@0*@1", gU, r)'.format(scale))
            self.modelBuilder.factory_('expr::gU4("@0*@0*@0*@0*@1", gU, r)')
          else:
            self.modelBuilder.factory_('expr::gU2("-{}*@0*@0", gU)'.format(scale))
            self.modelBuilder.factory_('expr::gU4("@0*@0*@0*@0", gU)')
        else:
          if self.grid:
            self.modelBuilder.factory_('expr::gU2("-{}*(((@0/fabs(@0))*@0)**0.5)*@1", mu, r)'.format(scale))
            self.modelBuilder.factory_('expr::gU4("@0*@1", mu, r)')
          else:
            self.modelBuilder.factory_('expr::gU2("-{}*(((@0/fabs(@0))*@0)**0.5)", mu)'.format(scale))
            self.modelBuilder.factory_('expr::gU4("@0", mu)')

    def getYieldScale(self, bin_, process):

        scalings = []

        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' not in process: 
            scalings.append('gU4')
        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' in process:
            scalings.append('gU2')
        #if "ggH_t_1200" in process:
        #    scalings.append('ggH_t_scale')
        #if "ggH_b_1200" in process:
        #    scalings.append('ggH_b_scale')
        #if "ggH_i_1200" in process:
        #    scalings.append('ggH_i_scale')

        if scalings:
            scaling = '_'.join(scalings)

            print 'Scaling %s/%s as %s' % (bin_, process,scaling)
            return scaling
        else:
            return 1

VLQ = VLQ()
