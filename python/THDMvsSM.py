#!/usr/bin/env python
import os
from collections import defaultdict
import math

import numpy as np

import ROOT

from HiggsAnalysis.CombinedLimit.PhysicsModel import *


class THDMvsSMHiggsModel(PhysicsModel):

    _binning = {
            "BP1_Type2": {
                "mH": np.arange(125., 1075., 50.),
                "tanb": np.arange(1.5, 51.5, 1.),
                },
            "BP1_Type1": {
                "mH": np.arange(125., 1075., 50.),
                "tanb": np.arange(1.5, 51.5, 1.),
                },
            "FixedMass_Type2": {
                "cos_betal": np.arange(-0.51, 0.53, 0.02),
                "tanb": np.arange(1., 51., 2.),
                },
    }


    def __init__(self):
        super(THDMvsSMHiggsModel, self).__init__()
        # Add containers to store set of defined signal processes,
        # the associated systematic uncertainties and nuisance parameters
        self.PROC_SETS = []
        self.SYST_DICT = defaultdict(list)
        self.NUISANCES = set()
        self.filename = ''
        self.debug_output = None
        self.x_variable = "mH"

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith('filePrefix='):
                filePrefix = po.replace('filePrefix=', '')
                print 'Set file prefix to: %s' % filePrefix
            if po.startswith('modelFile='):
                modelFile = po.replace('modelFile=', '')
                print "Importing scenario from '%s''"%(modelFile)
            if po.startswith('MSSM-NLO-Workspace='):
                self.ggHatNLO = po.replace('MSSM-NLO-Workspace=', '')
                print "Using %s for MSSM ggH NLO reweighting"%self.ggHatNLO
            if po.startswith('debug-output='):
                debug_name = po.replace('debug-output=', '')
                self.debug_output = ROOT.TFile.Open(debug_name, "recreate")
                print "Using %s as debug output file"%debug_name
        self.filename = os.path.join(filePrefix, modelFile)
        if "FixedMass" in modelFile.replace(".root", ""):
            self.x_variable = "cos_betal"

    def setModelBuilder(self, modelBuilder):
        """Used to load quantities in empty workspace."""
        # First call the parent class implementation
        super(THDMvsSMHiggsModel, self).setModelBuilder(modelBuilder)
        # Function to implement the histograms of (mA, tanb) dependent quantities
        self.buildModel()

    def buildModel(self):
        if self.x_variable == "mH":
            mass = ROOT.RooRealVar(self.x_variable, 'm_{H} [GeV]', 160.)
        elif self.x_variable == "cos_betal":
            mass = ROOT.RooRealVar(self.x_variable, 'cos(#beta-#alpha)', 0.1)
        tanb = ROOT.RooRealVar('tanb', 'tan#beta', 5.5)
        pars = [mass, tanb]

        procs = ["h", "H", "A"]

        # Open input root file to read the histograms
        rf = ROOT.TFile(self.filename, "read")

        for X in procs:
            if "H" == X: # don't create histogram for 'H' in cases, where its mass is a model-parameter
                continue
            self.doHistFunc("m_{}".format(X), rf.Get("m_{}".format(X)), pars) # syntax: Higgs-Boson, mass attribute, parameters

        for X in procs:
            self.doHistFunc("br_{}tautau".format(X), rf.Get("br_{}tautau".format(X)), pars) # syntax: Higgs-Boson, xsec attribute, parameters, production mode

            self.doHistFunc("xs_gg{}".format(X), rf.Get("xs_gg{}".format(X)), pars) # syntax: Higgs-Boson, xsec attribute, parameters, production mode
            self.doHistFunc("xs_bb{}".format(X), rf.Get("xs_bb{}".format(X)), pars) # syntax: Higgs-Boson, xsec attribute, parameters, production mode

            # self.doHistFuncFromXsecTools(X, "yukawa_top", pars)
            # self.doHistFuncFromXsecTools(X, "yukawa_bottom", pars)
            for loopcontrib in ['t', 'b']:
                # self.doHistFuncForYukawas(X, loopcontrib, pars)
                self.doHistFunc("Y{}_MSSM_{}".format(loopcontrib, X), rf.Get("g{}_{}".format(loopcontrib, X)), pars)

            self.add_ggH_at_NLO('xs_gg{X}{LC}', X)  # TODO: Maybe reenable ggH pT reweighting at NLO

            # ggH scale uncertainty
            self.doAsymPowSystematic(X, "xs", pars, "gg", "scale")
            # # ggH pdf+alpha_s uncertainty
            # self.doAsymPowSystematic(X, "xsec", pars, "gg", "pdfas")
            # # bbH total uncertainty
            # self.doAsymPowSystematic(X, "xsec", pars, "bb", "total")

            for loopcontrib in ['t','b','i']:
                self.SYST_DICT['xs_gg%s_%s' % (X, loopcontrib)].append('systeff_xs_gg%s_2HDM_scale' %X)

            # Make a note of what we've built, will be used to create scaling expressions later
            self.PROC_SETS.append('bb%s'%X)
            self.PROC_SETS.extend(['gg%s_t'%X, 'gg%s_b'%X, 'gg%s_i'%X])

        # And the SM terms
        self.PROC_SETS.extend(['ggH125', 'qqH125', 'ZH125', 'WH125', 'bbH125'])
        if self.debug_output:
            self.debug_output.Close()
        return

    def preProcessNuisances(self,nuisances):
        doParams = set()
        for bin in self.DC.bins:
            for proc in self.DC.exp[bin].keys():
                if self.DC.isSignal[proc]:
                    scaling = 'scaling_%s' % proc
                    print(scaling)
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

        # MSSMvsSM
        self.modelBuilder.doVar("x[1,0,1]")
        self.modelBuilder.out.var('x').setConstant(True)
        self.modelBuilder.factory_("expr::not_x(\"(1-@0)\", x)")
        self.sigNorms = { True:'x', False:'not_x' }

        # Create variable to set processes to zero
        self.modelBuilder.doVar('z[0,0,1]')
        self.modelBuilder.out.var('z').setConstant(True)

        self.modelBuilder.doSet('POI', 'r') # TODO: Check what this does and if it is correct

        # We don't intend on actually floating these in any fits...
        self.modelBuilder.out.var(self.x_variable).setConstant(True) # either mA or mHp
        self.modelBuilder.out.var('tanb').setConstant(True)

        # bsm_proc_match = "(gg(A|H|h)_(t|i|b)|bb(A|H|h))"
        bsm_proc_match = "(gg(A|H|h)_(t|b|i)|bb(A|H|h))"
        # Go over set of defined signal processes and scale the processes with the
        # correct normalizations.
        for proc in self.PROC_SETS:
            # Set up list of terms used in the scaling of the considered process
            terms = []
            # Get Higgs boson type from process name
            X = proc.split('_')[0].replace('gg','').replace('bb','')
            # For standard modell processes add flag to be able to deactivate them
            if "H125" in proc:
                terms = [self.sigNorms[False]]
            # In case of BSM signals scale the process with cross section and BR from
            # model prediction, free parameter r and flag it as BSM signal
            elif re.match(bsm_proc_match, proc): # not SM-like BSMSCALAR: either h or H
                terms.extend([
                    'xs_%s' % proc, # %proc.rstrip("_t"),  # TODO: Temp. fix for using only top contribution
                    'br_%stautau'%X,
                    'r',
                    self.sigNorms[True],
                    ])
            else:
                terms = ["z"]
            print(terms)

            # Check if a term was added that is associated with
            # a systematic uncertainty. If this is the case add
            # the systematic uncertainty to the scaling of the process
            extra = []
            for term in terms:
                if term in self.SYST_DICT:
                    extra += self.SYST_DICT[term]
            terms.extend(extra)
            print(proc, terms)
            # Add scaling function for the process to the workspace
            print('prod::scaling_%s(%s)'%(proc,','.join(terms)))
            self.modelBuilder.factory_('prod::scaling_%s(%s)'%(proc,','.join(terms)))
            self.modelBuilder.out.function('scaling_%s'%proc).Print('')

    def doAsymPowSystematic(self, higgs, quantity, varlist, production, uncertainty):
        # create AsymPow rate scaler given two TH2 inputs corresponding to kappa_hi and kappa_lo
        name = "{q}_{pr}{h}".format(q=quantity, pr=production, h=higgs)
        param = name + "_2HDM_" + uncertainty
        self.modelBuilder.doVar('%s[0,-7,7]'%param)
        param_var = self.modelBuilder.out.var(param)
        systname = "systeff_%s"%param

        # Open input root file to read the histograms
        rf = ROOT.TFile(self.filename, "read")
        hist_hi = rf.Get("{name}_{unc}_up".format(name=name, unc=uncertainty))
        hist_lo = rf.Get("{name}_{unc}_down".format(name=name, unc=uncertainty))
        # Currently the histograms in stored in the root file only represent
        # the value of the uncertainty. For a meaningful result the
        # prediction for the quantity needs to be added.
        hist_nom = rf.Get(name)
        hist_hi.Add(hist_nom)
        hist_hi.Divide(hist_nom)
        hist_lo.Add(hist_nom)
        hist_lo.Divide(hist_nom)

        self.NUISANCES.add(param)
        hi = self.doHistFunc('%s_hi'%systname, hist_hi, varlist)
        lo = self.doHistFunc('%s_lo'%systname, hist_lo, varlist)
        asym = ROOT.AsymPow(systname, '', lo, hi, param_var)
        self.modelBuilder.out._import(asym)
        return self.modelBuilder.out.function(systname)

    def getYieldScale(self,bin,process):
        if self.DC.isSignal[process]:
            scaling = 'scaling_%s' % process
            print 'Scaling %s/%s as %s' % (bin, process, scaling)
            return scaling
        else:
            return 1

    def doHistFunc(self, name, hist, varlist):
        if self.debug_output:
            self.debug_output.cd()
            hist.Write()
        dh = ROOT.RooDataHist('dh_%s'%name, 'dh_%s'%name, ROOT.RooArgList(*varlist), ROOT.RooFit.Import(hist))
        hfunc = ROOT.RooHistFunc(name, name, ROOT.RooArgSet(*varlist), dh)
        self.modelBuilder.out._import(hfunc, ROOT.RooFit.RecycleConflictNodes())
        return self.modelBuilder.out.function(name)

    def add_ggH_at_NLO(self, name, X):
        """Add separate scaling terms for each loop contribution in gluon fusion.

        Currently use only the SM fractions for the full range. To be replaced
        proper rescaling once the Yukawa couplings are available.
        """  # TODO
        fractions_sm = ROOT.TFile.Open(self.ggHatNLO, 'read')
        w_sm = fractions_sm.Get("w")

        template_X = X
        for loopcontrib in ['t','b','i']:
            func = w_sm.function('gg{X}_{LC}_MSSM_frac'.format(X=template_X, LC=loopcontrib))
            getattr(self.modelBuilder.out, 'import')(func, ROOT.RooFit.RecycleConflictNodes())
            self.modelBuilder.out.factory('prod::%s(%s,%s)' % (name.format(X=X, LC="_"+loopcontrib), name.format(X=X, LC=""), "gg%s_%s_MSSM_frac" % (template_X,loopcontrib))) #multiply t,b,i fractions with xsec at NNLO


    def doHistFuncForQQH(self, varlist):
        """Function to add scaling for qqh process of SM like Higgs boson.

        This is not needed if high-mass categorisation is used as there is no sensitivity to
        the qqh process in this categorisation.
        """
        # TODO: Is this description correct?
        raise NotImplementedError

    def doHistFuncForYukawas(self, X, contrib, varlist, thdm_type2=True):
        """Calculate Yukawa coupling from known value of cos(beta-alpha) and tanb value in grid.

        Temporary solution until the yukawa couplings will be available in the histograms.
        """
        name = "Y{c}_MSSM_{h}".format(c=contrib, h=X)
        x_parname = varlist[0].GetName()
        x_binning = self._binning[os.path.basename(self.filename).replace(".root", "")][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self._binning[os.path.basename(self.filename).replace(".root", "")][y_parname]

        # TODO: name needs to be fixed
        y_step = (y_binning[1]-y_binning[0])/2.
        hist = ROOT.TH2D(name, name, len(x_binning)-1, x_binning, len(y_binning)-1, y_binning)
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):
                beta = math.atan(hist.GetYaxis().GetBinCenter(i_y+1))
                alpha = beta - math.acos(0.01) # cos(beta-alpha) = 0.01

                if X == "h":
                    if contrib == "t":
                        value = math.cos(alpha) / math.sin(beta)
                    elif contrib == "b":
                        value = - math.sin(alpha) / math.cos(beta) if thdm_type2 else math.cos(alpha) / math.sin(beta)
                    else:
                        raise RuntimeError("Given contribution '{}' not known.".format(contrib))
                elif X == "H":
                    if contrib == "t":
                        value = math.sin(alpha) / math.sin(beta)
                    elif contrib == "b":
                        value = math.cos(alpha) / math.cos(beta) if thdm_type2 else math.sin(alpha) / math.sin(beta)
                    else:
                        raise RuntimeError("Given contribution '{}' not known.".format(contrib))
                elif X == "A":
                    if contrib == "t":
                        value = - 1. / math.tan(beta) # scaling as -cot(beta)
                    elif contrib == "b":
                        value = - math.tan(beta) if thdm_type2 else 1. / math.tan(beta)
                    else:
                        raise RuntimeError("Given contribution '{}' not known.".format(contrib))
                else:
                    raise RuntimeError("Given higgs boson '{}' not knwon.".format(X))

                hist.SetBinContent(i_x+1, i_y+1, value)

        return self.doHistFunc(name, hist, varlist)


THDMvsSM = THDMvsSMHiggsModel()
