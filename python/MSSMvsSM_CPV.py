from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from CombineHarvester.MSSMvsSMRun2Legacy.mssm_xs_tools import mssm_xs_tools
import CombineHarvester.CombineTools.plotting as plot

import os
import ROOT
import sys
import pprint
import numpy as np
import itertools
import re
from collections import defaultdict
from array import array

# tmp trick
CPV_to_classic = {
    'H1' : 'h',
    'H2' : 'H',
    'H3' : 'A',
    'h' : 'h',
    'H' : 'H',
    'A' : 'A',
}

class MSSMvsSMHiggsModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        ROOT.gROOT.SetBatch(ROOT.kTRUE)
        plot.ModTDRStyle(l=0.13, b=0.10, r=0.19)
        ROOT.gStyle.SetNdivisions(510, "Z")
        self.filePrefix = ''
        self.modelFile = ''
        self.filename = ''
        self.scenario = ''
        self.energy = ''
        self.ggHatNLO = ''
        self.mssm_inputs = None
        self.minTemplateMass = None
        self.maxTemplateMass = None
        self.quantity_map = {
            "mass"          : {"name" : "m{HIGGS}", "access" : "{HIGGS}"},
            "br"            : {"name" : "br_{HIGGS}tautau", "access": "{HIGGS}->tautau", "access2" : "br_{HIGGS}_tautau"},
            "xsec"          : {"name" : "xs_{PROD}{HIGGS}", "access" : "{PROD}->{HIGGS}", "access2" : "xs_{PROD}_{HIGGS}"},
            # In contrast to the CP conserving scenarios we do not provide relative Yukawa couplings as those do have a real and imaginary contribution to both the vector and axial-vector component. On the other hand, we do provide interference factors as in such channels large destructive interferences between H2 and H3 appear.
            "yukawa_top"    : {"name" : "Yt_MSSM_{HIGGS}", "access" : "rescale_gt_{HIGGS}"},
            "yukawa_bottom" : {"name" : "Yb_MSSM_{HIGGS}", "access" : "rescale_gb_{HIGGS}"},
            "yukawa_deltab" : {"name" : "Ydeltab_MSSM", "access" : "rescale_deltab"},
            "int_bb_tautau_H1" : {"name" : "int_bb_tautau_{HIGGS}", "access" : "int_bb_tautau_{HIGGS}"},
            "int_bb_tautau_H2" : {"name" : "int_bb_tautau_{HIGGS}", "access" : "int_bb_tautau_{HIGGS}"},
            "int_bb_tautau_H3" : {"name" : "int_bb_tautau_{HIGGS}", "access" : "int_bb_tautau_{HIGGS}"},
            "int_gg_tautau_H1" : {"name" : "int_gg_tautau_{HIGGS}", "access" : "int_gg_tautau_{HIGGS}"},
            "int_gg_tautau_H2" : {"name" : "int_gg_tautau_{HIGGS}", "access" : "int_gg_tautau_{HIGGS}"},
            "int_gg_tautau_H3" : {"name" : "int_gg_tautau_{HIGGS}", "access" : "int_gg_tautau_{HIGGS}"},
        }
        self.uncertainty_map = {
            "ggscale" : "::scale{VAR}",
            "ggpdfas" : "::pdfas{VAR}",
            "bbtotal" : "::{VAR}",
        }
        self.binning =  {
            "mh1125_CPV": {
                "tanb" : np.arange(1.0, 20.0, 1.0),
                "mHp" :   np.arange(130.0, 1500.0, 5.0),
            },
        }
        self.PROC_SETS = []
        self.SYST_DICT = defaultdict(list)
        self.NUISANCES = set()
        self.scaleforh = 1.0
        self.is_CPV = True # False by default later, to check how to set it

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
                self.scenario = self.modelFile.replace('.root','').replace('_%s'%self.energy,'').replace('_fixed','').replace('_original','').replace('_intermediate','')
                print "Importing scenario '%s' for sqrt(s) = '%s TeV' and '%s' data-taking period from '%s'"%(self.scenario, self.energy, self.era, self.modelFile)

            if po.startswith('MSSM-NLO-Workspace='):
                self.ggHatNLO = po.replace('MSSM-NLO-Workspace=', '')
                print "Using %s for MSSM ggH NLO reweighting"%self.ggHatNLO

            if po.startswith('minTemplateMass='):
                self.minTemplateMass = float(po.replace('minTemplateMass=', ''))
                print "Lower limit for mass histograms: {MINMASS}".format(MINMASS=self.minTemplateMass)

            if po.startswith('maxTemplateMass='):
                self.maxTemplateMass = float(po.replace('maxTemplateMass=', ''))
                print "Upper limit for mass histograms: {MAXMASS}".format(MAXMASS=self.maxTemplateMass)

            if po.startswith('scaleforh='):
                self.scaleforh = float(po.replace('scaleforh=',''))
                print "Additional scale for the light scalar h: {SCALE}".format(SCALE=self.scaleforh)

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
                value = getattr(self.mssm_inputs, quantity)(accesskey, x, y)
                if quantity == 'mass' and self.minTemplateMass:
                    if value < self.minTemplateMass:
                        print "[WARNING]: Found a value for {MH} below lower mass limit: {VALUE} < {MINMASS} for {XNAME} = {XVALUE}, {YNAME} = {YVALUE}. Setting it to limit".format(
                            MH=name,
                            VALUE=value,
                            MINMASS=self.minTemplateMass,
                            XNAME=x_parname,
                            XVALUE=x,
                            YNAME=y_parname,
                            YVALUE=y)
                        value = self.minTemplateMass
                if quantity == 'mass' and self.maxTemplateMass:
                    if value > self.maxTemplateMass:
                        print "[WARNING]: Found a value for {MH} above upper mass limit: {VALUE} > {MINMASS} for {XNAME} = {XVALUE}, {YNAME} = {YVALUE}. Setting it to limit".format(
                            MH=name,
                            VALUE=value,
                            MINMASS=self.maxTemplateMass,
                            XNAME=x_parname,
                            XVALUE=x,
                            YNAME=y_parname,
                            YVALUE=y)
                        value = self.maxTemplateMass
                hist.SetBinContent(i_x+1, i_y+1, value)
        return self.doHistFunc(name, hist, varlist)

    def doHistFuncForQQH(self, varlist):
        # Computing scaling function for qqh contribution (little Higgs) in context of MSSM
        name  = "qqH1_MSSM"
        if not self.is_CPV:
            accesskey = self.quantity_map['yukawa_top']['access'].format(HIGGS='H2')
        accesskey_br = self.quantity_map['br']['access2'].format(HIGGS='H1')
        print "Computing 'qqH1' scaling function from model file..."

        x_parname = varlist[0].GetName()
        x_binning = self.binning[self.scenario][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self.binning[self.scenario][y_parname]

        F = ROOT.TFile.Open(self.filename, "read")
        if not self.is_CPV:
            g_Htt_hist = F.Get(accesskey)
        br_htautau_hist = F.Get(accesskey_br)
        br_htautau_SM_125 = 0.06272 # Value for 125 GeV SM Higgs from YR4

        hist = ROOT.TH2D(name, name, len(x_binning)-1, x_binning, len(y_binning)-1, y_binning)
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):
                beta = np.arctan(y)
                if not self.is_CPV:
                    g_Htt = g_Htt_hist.Interpolate(x,y)
                br_htautau = br_htautau_hist.Interpolate(x,y)
                if not self.is_CPV:
                    sin_alpha = g_Htt * np.sin(beta)
                    if abs(sin_alpha) > 1:
                        sin_alpha  = np.sign(sin_alpha)
                    alpha = np.arcsin(sin_alpha)
                if not self.is_CPV:
                    value = np.sin(beta-alpha)**2 * br_htautau / br_htautau_SM_125 # (g_HVV)**2 * br_htautau / br_htautau_SM_125
                else:
                    value = br_htautau / br_htautau_SM_125
                value *= self.scaleforh # additional manual rescaling of light scalar h (default is 1.0)
                hist.SetBinContent(i_x+1, i_y+1, value)
        print "\trescale values range from",hist.GetMinimum(),"to",hist.GetMaximum()
        canv = ROOT.TCanvas()
        canv.cd()
        hist.SetContour(2000)
        hist.GetXaxis().SetTitle('m_{H^#pm}')
        hist.GetYaxis().SetTitle('tan#beta')
        hist.Draw("colz")
        canv.SetLogx()
        canv.Update()
        canv.SaveAs("qqH1_MSSM_%s.pdf"%self.scenario)
        canv.SaveAs("qqH1_MSSM_%s.png"%self.scenario)

        return self.doHistFunc(name, hist, varlist)

    def doHistFuncForGGH(self, varlist):
        # Computing scaling function for ggh contribution (little Higgs) in context of MSSM
        name  = "ggH1_MSSM"
        accesskey_xs = self.quantity_map['xsec']['access2'].format(HIGGS='H1',PROD='gg')
        accesskey_br = self.quantity_map['br']['access2'].format(HIGGS='H1')
        print "Computing 'ggH1' scaling function from model file..."

        x_parname = varlist[0].GetName()
        x_binning = self.binning[self.scenario][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self.binning[self.scenario][y_parname]

        F = ROOT.TFile.Open(self.filename, "read")
        xs_ggh_hist = F.Get(accesskey_xs)
        br_htautau_hist = F.Get(accesskey_br)
        br_htautau_SM_125 = 0.06272 # Value for 125 GeV SM Higgs from YR4
        xs_ggh_SM_125 = 48.58 #Value for 125 GeV SM Higgs from YR4

        hist = ROOT.TH2D(name, name, len(x_binning)-1, x_binning, len(y_binning)-1, y_binning)
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):
                xs_ggh = xs_ggh_hist.Interpolate(x,y)
                br_htautau = br_htautau_hist.Interpolate(x,y)
                value =  xs_ggh / xs_ggh_SM_125 * br_htautau / br_htautau_SM_125 # xs * BR / (xs * BR of SM 125)
                value *= self.scaleforh # additional manual rescaling of light scalar h (default is 1.0)
                hist.SetBinContent(i_x+1, i_y+1, value)
        print "\trescale values range from",hist.GetMinimum(),"to",hist.GetMaximum()
        canv = ROOT.TCanvas()
        canv.cd()
        hist.SetContour(2000)
        hist.GetXaxis().SetTitle('m_{H^#pm}')
        hist.GetYaxis().SetTitle('tan#beta')
        hist.Draw("colz")
        canv.SetLogx()
        canv.Update()
        canv.SaveAs("ggH1_MSSM_%s.pdf"%self.scenario)
        canv.SaveAs("ggH1_MSSM_%s.png"%self.scenario)

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

    def add_ggH_at_NLO_CPV(self, X):
        if not self.is_CPV:
            return False

        print("[INFO] Adding aditional terms for mssm ggh NLO reweighting for {}.".format(X))

        fractions_sm = ROOT.TFile(self.ggHatNLO, 'READ')
        w_sm = fractions_sm.Get("w")
        
        to_import = []
        
        mH = w_sm.var("mh")
        mH.SetName("m{}".format(X))
        to_import.append(mH)
        
        for LC in ["t", "b", "i"]:
            for _type in ["SM_frac", "SM_xsec", "2HDM_xsec"]:
                func = w_sm.function("ggh_{}_{}".format(LC, _type))
                # if _type == "SM_xsec":
                #     func.SetName("xs_gg{}_{}".format(X, LC))
                # else:
                #     func.SetName("gg{}_{}_{}".format(X, LC, _type).replace("_SM_frac", "_frac"))
                func.SetName("gg{}_{}_{}".format(X, LC, _type).replace("_SM_frac", "_frac"))
                print("Imported {} as {}".format("ggh_{}_{}".format(LC, _type), "gg{}_{}_{}".format(X, LC, _type).replace("_SM_frac", "_frac")))
                to_import.append(func)
        
        func = w_sm.function("ggh_SM_xsec")
        func.SetName("gg{}_SM_xsec".format(X))
        print("Imported {} as {}".format("ggh_SM_xsec", "gg{}_SM_xsec".format(X)))
        to_import.append(func)
        
        for i in to_import:
            getattr(self.modelBuilder.out, 'import')(i, ROOT.RooFit.RecycleConflictNodes())

        fractions_sm.Close()

        for loopcontrib in ["t", "b", "i"]:
            print('prod::%s(%s,%s)' % ('xs_gg{X}{LC}'.format(X=X, LC="_"+loopcontrib), 'xs_gg{X}{LC}'.format(X=X, LC=""), "gg%s_%s_frac" % (X,loopcontrib)))
            try:
                self.modelBuilder.out.factory('prod::%s(%s,%s)' % ('xs_gg{X}{LC}'.format(X=X, LC="_"+loopcontrib), 'xs_gg{X}{LC}'.format(X=X, LC=""), "gg%s_%s_frac" % (X,loopcontrib)))
            except:
                print('prod::%s(%s,%s)' % ('xs_gg{X}{LC}'.format(X=X, LC="_"+loopcontrib), 'xs_gg{X}{LC}'.format(X=X, LC=""), "gg%s_%s_frac" % (X,loopcontrib)))
                import pdb; pdb.set_trace()
        #'xs_gg{X}{LC}' = name.format(X=X, LC="") x "gg%s_%s_MSSM_frac" % (X,loopcontrib)
        return True

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
        self.modelBuilder.out.var('mHp').setConstant(True)
        self.modelBuilder.out.var('tanb').setConstant(True)

        for proc in self.PROC_SETS:
            X = None
            if re.match("(gg(H3|H2|H1)_(t|i|b)|bb(H3|H2|H1))", proc):
                X = proc.split('_')[0].replace('gg','').replace('bb','')
                terms = ['xs_%s' %proc, 'br_%stautau'%X]
                terms += ['r']
                terms += [self.sigNorms[True]]
            elif proc == 'qqH1':
                terms = [self.sigNorms[True], 'r', 'qqH1_MSSM']
            elif proc == 'ggH1':
                terms = [self.sigNorms[True], 'r', 'ggH1_MSSM']
            else:
                terms = [self.sigNorms[False]]
            if self.is_CPV and X in ['H2', 'H3']:
                for xx in ['bb', 'gg']:
                    if xx in proc:
                        terms.append('expr::interference_{PROD}_{HIGGS}(\"1.0 + @0\", int_{PROD}_tautau_{HIGGS})'.format(PROD=xx, HIGGS=X))
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
        mHp = ROOT.RooRealVar('mHp', 'm_{H^#pm} [GeV]', 130.)
        tanb = ROOT.RooRealVar('tanb', 'tan#beta', 10.)
        pars = [mHp, tanb]

        self.mssm_inputs = mssm_xs_tools(self.filename, True, 1) # syntax: model filename, Flag for interpolation ('True' or 'False'), verbosity level

        self.doHistFuncForQQH(pars)
        self.PROC_SETS.append('qqH1')

        self.doHistFuncForGGH(pars)
        self.PROC_SETS.append('ggH1')

        for X in ['H1', 'H2', 'H3']:
            self.doHistFuncFromXsecTools(X, "mass", pars) # syntax: Higgs-Boson, mass attribute, parameters

        self.doHistFuncFromModelFile(None, "yukawa_deltab", pars)

        for X in ['H1', 'H2', 'H3']:
            for xx in ['bb', 'gg']:
                self.doHistFuncFromModelFile(X, 'int_{}_tautau_{}'.format(xx, X), pars)

            self.doHistFuncFromXsecTools(X, "br", pars) # syntax: Higgs-Boson, xsec attribute, parameters, production mode

            self.doHistFuncFromXsecTools(X, "xsec", pars, production="gg") # syntax: Higgs-Boson, xsec attribute, parameters, production mode
            self.doHistFuncFromXsecTools(X, "xsec", pars, production="bb") # syntax: Higgs-Boson, xsec attribute, parameters, production mode
            if self.is_CPV:
                self.add_ggH_at_NLO_CPV(X)
            else:
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


MSSMvsSM_CPV = MSSMvsSMHiggsModel()
