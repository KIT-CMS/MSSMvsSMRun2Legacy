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
        self.vlq_split = False

    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("no_interference"):
                self.interference = False
            if po.startswith("grid"):
                self.grid = True
            if po.startswith("mu"):
                self.mu = True
            if po.startswith("vlq_split"):
                self.vlq_split = True


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

        if self.vlq_split:
          self.modelBuilder.doVar('betaLstau[0,0,1]')
          poiNames.append('betaLstau')


        if self.mu:
          self.modelBuilder.doVar('mu[0,0,10]')
          poiNames.append('mu')

        self.modelBuilder.doSet('POI', ','.join(poiNames))

        if self.vlq_split:
          if not self.mu:
            if self.grid:
              self.modelBuilder.factory_('expr::bbgU2("-{}*@0*@0*@1", gU, r)'.format(scale))
              self.modelBuilder.factory_('expr::bbgU4("@0*@0*@0*@0*@1", gU, r)')
              self.modelBuilder.factory_('expr::bsgU4("@0*@0*@0*@0*@1*(@2/0.19)*(@2/0.19)", gU, r, betaLstau)')
              self.modelBuilder.factory_('expr::ssgU2("-{}*@0*@0*@1*(@2/0.19)*(@2/0.19)*(@2/0.19)*(@2/0.19)", gU, r, betaLstau)'.format(scale))
              self.modelBuilder.factory_('expr::ssgU4("@0*@0*@0*@0*@1*(@2/0.19)*(@2/0.19)", gU, r, betaLstau)')

            else:
              self.modelBuilder.factory_('expr::bbgU2("-{}*@0*@0", gU)'.format(scale))
              self.modelBuilder.factory_('expr::bbgU4("@0*@0*@0*@0", gU)')
              self.modelBuilder.factory_('expr::bsgU4("@0*@0*@0*@0*(@1/0.19)*(@1/0.19)", gU, betaLstau)')
              self.modelBuilder.factory_('expr::ssgU2("-{}*@0*@0*(@1/0.19)*(@1/0.19)*(@1/0.19)*(@1/0.19)", gU, betaLstau)'.format(scale))
              self.modelBuilder.factory_('expr::ssgU4("@0*@0*@0*@0*(@1/0.19)*(@1/0.19)", gU, betaLstau)')
          else:
            if self.grid:
              self.modelBuilder.factory_('expr::bbgU2("-{}*(((@0/fabs(@0))*@0)**0.5)*@1", mu, r)'.format(scale))
              self.modelBuilder.factory_('expr::bbgU4("@0*@1", mu, r)')
              self.modelBuilder.factory_('expr::bsgU4("@0*@1*(@2/0.19)*(@2/0.19)", mu, r, betaLstau)')
              self.modelBuilder.factory_('expr::ssgU2("-{}*(((@0/fabs(@0))*@0)**0.5)*@1*(@2/0.19)*(@2/0.19)*(@2/0.19)*(@2/0.19)", mu, r, betaLstau)'.format(scale))
              self.modelBuilder.factory_('expr::ssgU4("@0*@1*(@2/0.19)*(@2/0.19)", mu, r, betaLstau)')

            else:
              self.modelBuilder.factory_('expr::bbgU2("-{}*(((@0/fabs(@0))*@0)**0.5)", mu)'.format(scale))
              self.modelBuilder.factory_('expr::bbgU4("@0", mu)')
              self.modelBuilder.factory_('expr::bsgU4("@0*(@1/0.19)*(@1/0.19)", mu, betaLstau)')
              self.modelBuilder.factory_('expr::ssgU2("-{}*(((@0/fabs(@0))*@0)**0.5)*(@1/0.19)*(@1/0.19)*(@1/0.19)*(@1/0.19)", mu, betaLstau)'.format(scale))
              self.modelBuilder.factory_('expr::ssgU4("@0*(@1/0.19)*(@1/0.19)", mu, betaLstau)')



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

        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' not in process and not process.startswith("bb") and not process.startswith("bs") and not process.startswith("ss"): 
            scalings.append('gU4')
        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' in process and not process.startswith("bb") and not process.startswith("bs") and not process.startswith("ss"):
            scalings.append('gU2')
        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' not in process and process.startswith("bb"):
            scalings.append('bbgU4')
        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' in process and process.startswith("bb"):
            scalings.append('bbgU2')
        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' not in process and process.startswith("bs"):
            scalings.append('bsgU4')
        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' not in process and process.startswith("ss"):
            scalings.append('ssgU4')
        if ('VLQ_betaRd33_0' in process or 'VLQ_betaRd33_minus1' in process) and 'interference' in process and process.startswith("ss"):
            scalings.append('ssgU2')

        if scalings:
            scaling = '_'.join(scalings)

            print 'Scaling %s/%s as %s' % (bin_, process,scaling)
            return scaling
        else:
            return 1

VLQ = VLQ()
