from HiggsAnalysis.CombinedLimit.PhysicsModel import *

import os
import ROOT
import sys
import json
import pprint
import numpy as np
import itertools
import re
from collections import defaultdict
from array import array

class YtYbScan(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        self.XS_File = ''
        self.PROC_SETS = []
        self.NUISANCES = set()
        self.SYST_DICT = defaultdict(list)

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith('XS-Workspace='):
                self.XS_File = po.replace('XS-Workspace=', '')
                print "Using %s for XS inputs"%self.XS_File

    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self, modelBuilder)
        self.buildModel()

    def doAsymPowSystematic(self, higgs, production, uncertainty):
        name = 'xs_{PROD}{HIGGS}'.format(HIGGS=higgs, PROD=production)


        # create AsymPow rate scaler 
        param = name + '_' + uncertainty
        self.modelBuilder.doVar('%s[0,-7,7]'%param)
        param_var = self.modelBuilder.out.var(param)
        systname = "systeff_%s"%param

        print "Doing AsymPow systematic '%s'" %(param)

        self.NUISANCES.add(param)
        hi = self.modelBuilder.out.function('xs_%(production)s%(higgs)s_up' % vars()) 
        lo = self.modelBuilder.out.function('xs_%(production)s%(higgs)s_down' % vars()) 
        asym = ROOT.AsymPow(systname, '', lo, hi, param_var)
        self.modelBuilder.out._import(asym)
        return self.modelBuilder.out.function(systname)

    def add_XS(self, X):

        f_w_xs = ROOT.TFile.Open(self.XS_File, 'read')
        w_xs = f_w_xs.Get("w")

        func = w_xs.function('xs_gg%(X)s' % vars())
        func_up = w_xs.function('xs_gg%(X)s_up' % vars())
        func_down = w_xs.function('xs_gg%(X)s_down' % vars())
        getattr(self.modelBuilder.out, 'import')(func, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.modelBuilder.out, 'import')(func_up, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.modelBuilder.out, 'import')(func_down, ROOT.RooFit.RecycleConflictNodes())


        for loopcontrib in ['t','b','i']:
            func = w_xs.function('xs_gg%(X)s_%(loopcontrib)s' % vars())
            func_up = w_xs.function('xs_gg%(X)s_%(loopcontrib)s_up' % vars())
            func_down = w_xs.function('xs_gg%(X)s_%(loopcontrib)s_down' % vars())
            getattr(self.modelBuilder.out, 'import')(func, ROOT.RooFit.RecycleConflictNodes())
            getattr(self.modelBuilder.out, 'import')(func_up, ROOT.RooFit.RecycleConflictNodes())
            getattr(self.modelBuilder.out, 'import')(func_down, ROOT.RooFit.RecycleConflictNodes())

        func = w_xs.function('xs_bb%(X)s' % vars())
        func_up = w_xs.function('xs_bb%(X)s_up' % vars())
        func_down = w_xs.function('xs_bb%(X)s_down' % vars())
        getattr(self.modelBuilder.out, 'import')(func, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.modelBuilder.out, 'import')(func_up, ROOT.RooFit.RecycleConflictNodes())
        getattr(self.modelBuilder.out, 'import')(func_down, ROOT.RooFit.RecycleConflictNodes())
        
    def preProcessNuisances(self,nuisances):
        doParams = set()
        for bin in self.DC.bins:
            for proc in self.DC.exp[bin].keys():
                if self.DC.isSignal[proc]:
                    scaling = 'scaling_%s' % proc
                    params = self.modelBuilder.out.function(scaling).getParameters(ROOT.RooArgSet()).contentsString().split(',')
                    for param in params:
                        if param in self.NUISANCES:
                            doParams.add(param)
        for param in doParams:
            print 'Add nuisance parameter %s to datacard' % param
            nuisances.append((param,False, "param", [ "0", "1"], [] ) )

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("r[1]") # overall SF r -> can act as BR scaling if needed
        self.modelBuilder.doSet('POI', 'r')

        # We don't intend on actually floating these in any fits...
        self.modelBuilder.out.var('mA').setConstant(True)
        self.modelBuilder.out.var('mH').setConstant(True)

        bsm_proc_match = "(gg(A|H|h|H3|H2|H1)_(t|i|b)|bb(A|H|h|H3|H2|H1))(|_lowmass)"
 
        for proc in self.PROC_SETS:
            terms = []
            X = proc.split('_')[0].replace('gg','').replace('bb','')
            if re.match(bsm_proc_match, proc): 
                terms = ['xs_%s' %proc.replace("_lowmass","")]
                terms += ['r']

            # Now scan terms and add theory uncerts
            extra = []
            for term in terms:
                if term in self.SYST_DICT:
                    extra += self.SYST_DICT[term]
            terms += extra
            self.modelBuilder.factory_('prod::scaling_%s(%s)'%(proc,','.join(terms)))
            self.modelBuilder.out.function('scaling_%s'%proc).Print('')
            

    def getYieldScale(self,bin,process):
        if self.DC.isSignal[process]:
            scaling = 'scaling_%s' % process
            print 'Scaling %s/%s as %s' % (bin, process, scaling)
            return scaling
        else:
            return 1

    def buildModel(self):

        procs = ['H', 'A']

        for X in procs:

            mass = ROOT.RooRealVar('m%(X)s' % vars(), 'm_{%(X)s} [GeV]' % vars(), 100.)
            Yt = ROOT.RooRealVar('Yt_%(X)s' % vars(), 'Y_{t}^{%(X)s}' % vars(), 1.)
            Yb = ROOT.RooRealVar('Yb_%(X)s' % vars(), 'Y_{b}^{%(X)s}' % vars(), 1.)

            pars = [mass, Yt, Yb]

            # get cross sections from other workspace
            self.add_XS(X)

            # ggH total uncertainty
            self.doAsymPowSystematic(X, "gg", "total")

            for loopcontrib in ['t','b','i']:
                self.SYST_DICT['xs_gg%s_%s' % (X, loopcontrib)].append('systeff_xs_gg%s_total' %X)

            # bbH total uncertainty
            self.doAsymPowSystematic(X, "bb", "total")
            self.SYST_DICT['xs_bb%s' %X].append('systeff_xs_bb%s_total' %X)

            # Make a note of what we've built, will be used to create scaling expressions later
            self.PROC_SETS.append('bb%s'%X)
            self.PROC_SETS.extend(['gg%s_t'%X, 'gg%s_b'%X, 'gg%s_i'%X])

YtYbScan = YtYbScan()
