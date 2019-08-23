from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from CombineHarvester.MSSMvsSMRun2Legacy.mssm_xs_tools import mssm_xs_tools

import os
import ROOT
import sys
import pprint
import numpy as np
import itertools
import re
from collections import defaultdict

class MSSMvsSMHiggsModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
        self.filePrefix = ''
        self.modelFile = ''
        self.filename = ''
        self.scenario = ''
        self.energy = ''
        self.ggHatNLO = ''
        self.mssm_inputs = None
        self.quantity_map = {
            "mass"          : {"name" : "m{HIGGS}", "access" : "{HIGGS}"},
            "br"            : {"name" : "br_{HIGGS}tautau", "access": "{HIGGS}->tautau"},
            "xsec"          : {"name" : "xs_{PROD}{HIGGS}", "access" : "{PROD}->{HIGGS}"},
            "yukawa_top"    : {"name" : "Yt_MSSM_{HIGGS}", "access" : "rescale_gt_{HIGGS}"},
            "yukawa_bottom" : {"name" : "Yb_MSSM_{HIGGS}", "access" : "rescale_gb_{HIGGS}"},
            "yukawa_deltab" : {"name" : "Ydeltab_MSSM", "access" : "rescale_deltab"},
        }
        self.uncertainty_map = {
            "ggscale" : "::scale{VAR}",
            "ggpdfas" : "::pdfas{VAR}",
            "bbtotal" : "::{VAR}",
        }
        self.binning =  {
            "mh125": {
                "tanb" : np.arange(1.0, 61.0, 1.0),
                "mA" :   np.arange(70.0, 2605.0, 5.0),
            },
            "mh125_lc": {
                "tanb" : np.arange(1.0, 61.0, 1.0),
                "mA" :   np.arange(70.0, 2605.0, 5.0),
            },
            "mh125_ls": {
                "tanb" : np.arange(1.0, 61.0, 1.0),
                "mA" :   np.arange(70.0, 2605.0, 5.0),
            },
            "mh125_alignment": {
                "tanb" : np.arange(1.0, 20.25, 0.25),
                "mA" :   np.arange(120.0, 1005.0, 5.0),
            },
        }
        self.PROC_SETS = []
        self.SYST_DICT = defaultdict(list)
        self.NUISANCES = set()

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith('filePrefix='):
                self.filePrefix = po.replace('filePrefix=', '')
                print 'Set file prefix to: %s' % self.filePrefix

            if po.startswith('modelFile='):
                cfg= po.replace('modelFile=', '')
                cfgSplit = cfg.split(',')
                if len(cfgSplit) != 3:
                    raise RuntimeError, 'Model file argument %s should be in the format ENERGY,ERA,FILE' % cfg
                self.energy = cfgSplit[0]
                self.era = cfgSplit[1]
                self.modelFile = cfgSplit[2]
                self.scenario = self.modelFile.replace('.root','').replace('_%s'%self.energy,'')
                print "Importing scenario '%s' for sqrt(s) = '%s TeV' and '%s' data-taking period from '%s'"%(self.scenario, self.energy, self.era, self.modelFile)

            if po.startswith('MSSM-NLO-Workspace='):
                self.ggHatNLO = po.replace('MSSM-NLO-Workspace=', '')
                print "Using %s for MSSM ggH NLO reweighting"%self.ggHatNLO

        self.filename = os.path.join(self.filePrefix, self.modelFile)

    def setModelBuilder(self, modelBuilder):
        # First call the parent class implementation
        PhysicsModel.setModelBuilder(self, modelBuilder)
        # Function to implement the histograms of (mA, tanb) dependent quantities
        self.buildModel()

    def doHistFunc(self, name, hist, varlist):
        dh = ROOT.RooDataHist('dh_%s'%name, 'dh_%s'%name, ROOT.RooArgList(*varlist), ROOT.RooFit.Import(hist))
        hfunc = ROOT.RooHistFunc(name, name, ROOT.RooArgSet(*varlist), dh)
        self.modelBuilder.out._import(hfunc, ROOT.RooFit.RecycleConflictNodes())
        return self.modelBuilder.out.function(name)

    def doHistFuncFromXsecTools(self, higgs, quantity, varlist, production=None):
        # Translator mssm_xs_tools -> TH1D -> RooDataHist
        name  = self.quantity_map[quantity]['name']
        accesskey = self.quantity_map[quantity]['access']
        if production and higgs:
            name = name.format(HIGGS=higgs, PROD=production)
            accesskey = accesskey.format(HIGGS=higgs, PROD=production)
        elif higgs:
            name = name.format(HIGGS=higgs)
            accesskey = accesskey.format(HIGGS=higgs)
        print "Doing histFunc '%s' with '%s' key for quantity '%s' from mssm_xs_tools..." %(name, accesskey, quantity)

        x_parname = varlist[0].GetName()
        x_binning = self.binning[self.scenario][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self.binning[self.scenario][y_parname]

        hist = ROOT.TH2D(name, name, len(x_binning)-1, x_binning, len(y_binning)-1, y_binning) 
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):
                hist.SetBinContent(i_x+1, i_y+1, getattr(self.mssm_inputs, quantity)(accesskey, x, y))
        return self.doHistFunc(name, hist, varlist)

    def doHistFuncFromModelFile(self, higgs, quantity, varlist):
        name  = self.quantity_map[quantity]['name']
        accesskey = self.quantity_map[quantity]['access']
        if higgs:
            name = name.format(HIGGS=higgs)
            accesskey = accesskey.format(HIGGS=higgs)
        F = ROOT.TFile.Open(self.filename, "read")
        print "Doing histFunc '%s' with '%s' key for quantity '%s' from model file..." %(name, accesskey, quantity)
        hist = F.Get(accesskey)
        hist.SetName(name)
        return self.doHistFunc(name, hist, varlist)

    def doAsymPowSystematic(self, higgs, quantity, varlist, production, uncertainty):
        # Translator mssm_xs_tools -> TH1D -> RooDataHist -> Systematic
        name  = self.quantity_map[quantity]['name'].format(HIGGS=higgs, PROD=production)
        accesskey = self.quantity_map[quantity]['access'].format(HIGGS=higgs, PROD=production)
        uncertaintykey = self.uncertainty_map[production+uncertainty]

        # create AsymPow rate scaler given two TH2 inputs corresponding to kappa_hi and kappa_lo
        param = name + "_MSSM_" + uncertainty
        self.modelBuilder.doVar('%s[0,-7,7]'%param)
        param_var = self.modelBuilder.out.var(param)
        systname = "systeff_%s"%param

        x_parname = varlist[0].GetName()
        x_binning = self.binning[self.scenario][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self.binning[self.scenario][y_parname]

        hist_hi = ROOT.TH2D(systname+"_hi", systname+"_hi", len(x_binning)-1, x_binning, len(y_binning)-1, y_binning) 
        hist_lo = ROOT.TH2D(systname+"_lo", systname+"_lo", len(x_binning)-1, x_binning, len(y_binning)-1, y_binning) 
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):
                nominal  = getattr(self.mssm_inputs, quantity)(accesskey, x, y)
                value_hi = getattr(self.mssm_inputs, quantity)(accesskey+uncertaintykey.format(VAR='up'), x, y) 
                value_lo = getattr(self.mssm_inputs, quantity)(accesskey+uncertaintykey.format(VAR='down'), x, y) 
                if nominal == 0:
                    hist_hi.SetBinContent(i_x+1, i_y+1, 1.0)
                    hist_lo.SetBinContent(i_x+1, i_y+1, 1.0)
                else:
                    hist_hi.SetBinContent(i_x+1, i_y+1, (nominal+value_hi)/nominal)
                    hist_lo.SetBinContent(i_x+1, i_y+1, (nominal+value_lo)/nominal)
        print "Doing AsymPow systematic '%s' with '%s' key for quantity '%s' from mssm_xs_tools..." %(param, accesskey+uncertaintykey.format(VAR='up/down'), quantity)

        self.NUISANCES.add(param)
        hi = self.doHistFunc('%s_hi'%systname, hist_hi, varlist)
        lo = self.doHistFunc('%s_lo'%systname, hist_lo, varlist)
        asym = ROOT.AsymPow(systname, '', lo, hi, param_var)
        self.modelBuilder.out._import(asym)
        return self.modelBuilder.out.function(systname)

    def add_ggH_at_NLO(self, name, X):
        importstring = os.path.expandvars(self.ggHatNLO)+":w:gg{X}_{LC}_MSSM_frac" #import t,b,i fraction of xsec at NLO
        for loopcontrib in ['t','b','i']:
            getattr(self.modelBuilder.out, 'import')(importstring.format(X=X, LC=loopcontrib), ROOT.RooFit.RecycleConflictNodes())
            self.modelBuilder.out.factory('prod::%s(%s,%s)' % (name.format(X=X, LC="_"+loopcontrib), name.format(X=X, LC=""), "gg%s_%s_MSSM_frac" % (X,loopcontrib))) #multiply t,b,i fractions with xsec at NNLO

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
        self.modelBuilder.doVar("r[1,0,20]")

        #MSSMvsSM
        self.modelBuilder.doVar("x[1,0,1]")
        self.modelBuilder.out.var('x').setConstant(True)
        self.modelBuilder.factory_("expr::not_x(\"(1-@0)\", x)")
        self.sigNorms = { True:'x', False:'not_x' }

        self.modelBuilder.doSet('POI', 'r')
        
        # We don't intend on actually floating these in any fits...
        self.modelBuilder.out.var('mA').setConstant(True)
        self.modelBuilder.out.var('tanb').setConstant(True)

        for proc in self.PROC_SETS:
            if re.match("(gg(A|H|h)_(t|i|b)|bb(A|H|h))", proc):
                terms = ['xs_%s' %proc, 'br_%stautau'%proc.split('_')[0].replace('gg','').replace('bb','')]
                terms += ['r']
                terms += [self.sigNorms[True]]
            else:
                terms = [self.sigNorms[False]]
            # Now scan terms and add theory uncerts
            extra = []
            for term in terms:
                if term in self.SYST_DICT:
                    extra += self.SYST_DICT[term]
            terms += extra
            self.modelBuilder.factory_('prod::scaling_%s(%s)'%(proc,','.join(terms)))
            self.modelBuilder.out.function('scaling_%s'%proc).Print('')


    def getHiggsProdDecMode(self, bin, process):
        """Return a triple of (production, decay, energy)"""
        P = ''
        D = ''
        if "_" in process:
            (P, D) = process.split("_")
        else:
            raise RuntimeError, 'Expected signal process %s to be of the form PROD_DECAY' % process
        E = None
        for era in self.ERAS:
            if era in bin:
                if E: raise RuntimeError, "Validation Error: bin string %s contains multiple known energies" % bin
                E = era
        if not E:
                raise RuntimeError, 'Did not find a valid energy in bin string %s' % bin
        return (P, D, E)

    def getYieldScale(self,bin,process):
        if self.DC.isSignal[process]:
            scaling = 'scaling_%s' % process
            print 'Scaling %s/%s as %s' % (bin, process, scaling)
            return scaling
        else:
            return 1

    def buildModel(self):
        mA = ROOT.RooRealVar('mA', 'm_{A} [GeV]', 120.)
        tanb = ROOT.RooRealVar('tanb', 'tan#beta', 20.)
        pars = [mA, tanb]

        self.mssm_inputs = mssm_xs_tools(self.filename, True, 1) # syntax: model filename, Flag for interpolation ('True' or 'False'), verbosity level

        for X in ['h', 'H']:
            self.doHistFuncFromXsecTools(X, "mass", pars) # syntax: Higgs-Boson, mass attribute, parameters

        self.doHistFuncFromModelFile(None, "yukawa_deltab", pars)

        for X in ['h', 'H', 'A']:
            self.doHistFuncFromModelFile(X, "yukawa_top", pars)
            self.doHistFuncFromModelFile(X, "yukawa_bottom", pars)

            self.doHistFuncFromXsecTools(X, "br", pars) # syntax: Higgs-Boson, xsec attribute, parameters, production mode

            self.doHistFuncFromXsecTools(X, "xsec", pars, production="gg") # syntax: Higgs-Boson, xsec attribute, parameters, production mode
            self.doHistFuncFromXsecTools(X, "xsec", pars, production="bb") # syntax: Higgs-Boson, xsec attribute, parameters, production mode
            self.add_ggH_at_NLO('xs_gg{X}{LC}', X)

            # ggH scale uncertainty
            self.doAsymPowSystematic(X, "xsec", pars, "gg", "scale")
            # ggH pdf+alpha_s uncertainty
            self.doAsymPowSystematic(X, "xsec", pars, "gg", "pdfas")
            # bbH total uncertainty
            self.doAsymPowSystematic(X, "xsec", pars, "bb", "total")

            self.SYST_DICT['xs_gg%s' %X].append('systeff_xs_gg%s_MSSM_scale' %X)
            self.SYST_DICT['xs_gg%s' %X].append('systeff_xs_gg%s_MSSM_pdfas' %X)
            for loopcontrib in ['t','b','i']:
                self.SYST_DICT['xs_gg%s_%s' % (X, loopcontrib)].append('systeff_xs_gg%s_MSSM_scale' %X)
                self.SYST_DICT['xs_gg%s_%s' % (X, loopcontrib)].append('systeff_xs_gg%s_MSSM_pdfas' %X)

            self.SYST_DICT['xs_bb%s' %X].append('systeff_xs_bb%s_MSSM_total' %X)

            # Make a note of what we've built, will be used to create scaling expressions later
            self.PROC_SETS.append('bb%s'%X)
            self.PROC_SETS.extend(['gg%s_t'%X, 'gg%s_b'%X, 'gg%s_i'%X])

        # And the SM terms
        self.PROC_SETS.extend(['ggH125', 'qqH125', 'ZH125', 'WH125', 'ttH125'])


MSSMvsSM = MSSMvsSMHiggsModel()
