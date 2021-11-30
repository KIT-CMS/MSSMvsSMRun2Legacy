from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel
from HiggsAnalysis.CombinedLimit.SMHiggsBuilder import SMHiggsBuilder
import numpy as np

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
        
        if self.mu:
          self.modelBuilder.doVar('mu[0,0,10]')
          poiNames.append('mu')

        self.modelBuilder.doSet('POI', ','.join(poiNames))

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

        if scalings:
            scaling = '_'.join(scalings)

            print 'Scaling %s/%s as %s' % (bin_, process,scaling)
            return scaling
        else:
            return 1

VLQ = VLQ()
