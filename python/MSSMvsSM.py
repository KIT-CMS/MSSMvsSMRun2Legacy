from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from CombineHarvester.MSSMvsSMRun2Legacy.mssm_xs_tools import mssm_xs_tools

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

class MSSMvsSMHiggsModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
        self.filePrefix = ''
        self.modelFile = ''
        self.filename = ''
        self.scenario = ''
        self.energy = ''
        self.ggHatNLO = ''
        self.mssm_inputs = None
        self.sm_predictions = None
        self.debug_output = None
        self.minTemplateMass = None
        self.maxTemplateMass = None
        self.quantity_map = {
            "mass"             : {"method" :     "mass", "name" : "m{HIGGS}",                 "access" : "{HIGGS}"},
            "width"            : {"method" :    "width", "name" : "w{HIGGS}",                 "access" : "{HIGGS}"},
            "width_SM"         : {"method" :    "width", "name" : "w{HIGGS}_SM",              "access" : "HSM"},
            "br"               : {"method" :       "br", "name" : "br_{HIGGS}tautau",         "access" : "{HIGGS}->tautau"},
            "br_SM"            : {"method" :       "br", "name" : "br_{HIGGS}tautau_SM",      "access" : "HSM->tautau"},
            "xsec"             : {"method" :     "xsec", "name" : "xs_{PROD}{HIGGS}",         "access" : "{PROD}->{HIGGS}"},
            "xsec_SM"          : {"method" :     "xsec", "name" : "xs_{PROD}{HIGGS}_SM",      "access" : "{PROD}->{ADD}HSM"},
            "interference"     : {"method" :     "xsec", "name" : "int_{PROD}{HIGGS}_tautau", "access" : "int_{PROD}_tautau_{HIGGS}"},
            "yukawa_top"       : {"method" : "coupling", "name" : "Yt_MSSM_{HIGGS}",          "access" : "gt_{HIGGS}"},
            "yukawa_bottom"    : {"method" : "coupling", "name" : "Yb_MSSM_{HIGGS}",          "access" : "gb_{HIGGS}"},
            "yukawa_deltab"    : {"method" : "coupling", "name" : "Ydeltab_MSSM",             "access" : "deltab"},
            "yukawa_im_deltab" : {"method" : "coupling", "name" : "Yimdeltab_MSSM",           "access" : "im_deltab"},
        }
        self.uncertainty_map = {
            "ggscale" : "::scale{VAR}",
            "ggpdfas" : "::pdfas{VAR}",
            "bbtotal" : "::{VAR}",
        }
        self.binning =  {
            "hMSSM": {
                "tanb" : np.concatenate((np.arange(0.5, 6.0, 0.1), np.arange(6.0, 61.0, 1.0))),
                "mA" :   np.arange(130.0, 2605.0, 5.0),
            },
            "mh125": {
                "tanb" : np.concatenate((np.arange(0.5, 1.0, 0.1), np.arange(1.0, 10.0, 0.5), np.arange(10.0, 61.0, 1.0))),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0), np.arange(200.0, 320.0, 5.0), np.arange(320.0, 370, 1.0), np.arange(370.0, 2605.0, 5.0))),
            },
            "mh125EFT": {
                "tanb" : np.arange(1.0, 10.25, 0.25),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0),np.arange(200.0, 320.0, 5.0),np.arange(320.0, 370.0, 1.0), np.arange(370.0, 3005.0, 5.0))),
            },
            "mh125_lc": {
                "tanb" : np.concatenate((np.arange(0.5, 1.0, 0.1), np.arange(1.0, 10.0, 0.5), np.arange(10.0, 61.0, 1.0))),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0), np.arange(200.0, 320.0, 5.0), np.arange(320.0, 370, 1.0), np.arange(370.0, 2605.0, 5.0))),
            },
            "mh125EFT_lc": {
                "tanb" : np.arange(1.0, 10.25, 0.25),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0),np.arange(200.0, 320.0, 5.0),np.arange(320.0, 370.0, 1.0), np.arange(370.0, 3005.0, 5.0))),
            },
            "mh125_ls": {
                "tanb" : np.concatenate((np.arange(0.5, 1.0, 0.1), np.arange(1.0, 10.0, 0.5), np.arange(10.0, 61.0, 1.0))),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0), np.arange(200.0, 320.0, 5.0), np.arange(320.0, 370, 1.0), np.arange(370.0, 2605.0, 5.0))),
            },
            "mh125_align": {
                "tanb" : np.arange(1.0, 20.25, 0.25),
                "mA" :   np.concatenate((np.arange(120.0, 600.0, 1.0), np.arange(600.0, 1005.0, 5.0))),
            },
            "mh125_muneg_1": {
                "tanb" : np.concatenate((np.arange(0.5, 1.0, 0.1), np.arange(1.0, 10.0, 0.5), np.arange(10.0, 61.0, 1.0))),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0), np.arange(200.0, 320.0, 5.0), np.arange(320.0, 370, 1.0), np.arange(370.0, 2605.0, 5.0))),
            },
            "mh125_muneg_2": {
                "tanb" : np.concatenate((np.arange(0.5, 1.0, 0.1), np.arange(1.0, 10.0, 0.5), np.arange(10.0, 61.0, 1.0))),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0), np.arange(200.0, 320.0, 5.0), np.arange(320.0, 370, 1.0), np.arange(370.0, 2605.0, 5.0))),
            },
            "mh125_muneg_3": {
                "tanb" : np.concatenate((np.arange(0.5, 1.0, 0.1), np.arange(1.0, 10.0, 0.5), np.arange(10.0, 61.0, 1.0))),
                "mA" :   np.concatenate((np.arange(70.0, 200.0, 1.0), np.arange(200.0, 320.0, 5.0), np.arange(320.0, 370, 1.0), np.arange(370.0, 2605.0, 5.0))),
            },
            "mHH125": {
                "tanb" : np.arange(5.0, 6.01, 0.01),
                "mHp" :   np.arange(150.0, 200.2, 0.2),
            },
            "mh1125_CPV": {
                "tanb" : np.arange(1.0, 20.25, 0.25),
                "mHp" :   np.concatenate((np.arange(130.0, 200.0, 1.0),np.arange(200.0, 320.0, 5.0),np.arange(320.0, 370.0, 1.0), np.arange(370.0, 1505.0, 5.0))),
            },
        }
        self.PROC_SETS = []
        self.SYST_DICT = defaultdict(list)
        self.NUISANCES = set()
        self.scaleforh = 1.0
        self.bsmscalar = ""
        self.smlike = "h"
        self.massparameter = "mA"
        self.replace_with_sm125 = True
        self.use_hSM_difference = False
        self.cpv_template_pairs = {"H1" : "h", "H2" : "H", "H3" : "A"}

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

                if self.scenario == "mHH125":
                    self.smlike = "H"
                    self.bsmscalar = "h"
                    self.massparameter = "mHp"
                elif self.scenario == "mh1125_CPV":
                    self.smlike = "H1"
                    self.bsmscalar = ""
                    self.massparameter = "mHp"
                else:
                    self.smlike = "h"
                    self.bsmscalar = "H"
                    self.massparameter = "mA"
                print "Chosen model-specific settings:"
                print "SM-like Higgs boson:",self.smlike
                print "BSM scalar Higgs boson:",self.bsmscalar
                print "Mass parameter in the plane:",self.massparameter

            if po.startswith('debug-output='):
                debug_name = po.replace('debug-output=', '')
                self.debug_output = ROOT.TFile.Open(debug_name, "recreate")
                print "Using %s as debug output file"%debug_name

            if po.startswith('MSSM-NLO-Workspace='):
                self.ggHatNLO = po.replace('MSSM-NLO-Workspace=', '')
                print "Using %s for MSSM ggH NLO reweighting"%self.ggHatNLO

            if po.startswith('replace-with-SM125='):
                self.replace_with_sm125 = bool(int(po.replace('replace-with-SM125=', ''))) # use either 1 or 0 for the choice
                print "Replacing with SM 125?",self.replace_with_sm125

            if po.startswith('sm-predictions='):
                sm_pred_path = po.replace('sm-predictions=','')
                self.sm_predictions = json.load(open(sm_pred_path,'r'))
                print "Using %s for SM predictions"%sm_pred_path

            if po.startswith('minTemplateMass='):
                self.minTemplateMass = float(po.replace('minTemplateMass=', ''))
                print "Lower limit for mass histograms: {MINMASS}".format(MINMASS=self.minTemplateMass)

            if po.startswith('maxTemplateMass='):
                self.maxTemplateMass = float(po.replace('maxTemplateMass=', ''))
                print "Upper limit for mass histograms: {MAXMASS}".format(MAXMASS=self.maxTemplateMass)

            if po.startswith('scaleforh='):
                self.scaleforh = float(po.replace('scaleforh=',''))
                print "Additional scale for the light scalar h: {SCALE}".format(SCALE=self.scaleforh)

            if po.startswith('hSM-treatment='):
                hSM_treatment = po.replace('hSM-treatment=', '')
                self.use_hSM_difference = hSM_treatment == "hSM-in-bg"
                print "Using (BSM - SM) difference for SM-like Higgs boson?",self.use_hSM_difference

        self.filename = os.path.join(self.filePrefix, self.modelFile)

    def setModelBuilder(self, modelBuilder):
        # First call the parent class implementation
        PhysicsModel.setModelBuilder(self, modelBuilder)
        # Function to implement the histograms of (mA, tanb) dependent quantities
        self.buildModel()

    def doHistFunc(self, name, hist, varlist):
        if self.debug_output:
            self.debug_output.cd()
            hist.Write()
        dh = ROOT.RooDataHist('dh_%s'%name, 'dh_%s'%name, ROOT.RooArgList(*varlist), ROOT.RooFit.Import(hist))
        hfunc = ROOT.RooHistFunc(name, name, ROOT.RooArgSet(*varlist), dh)
        self.modelBuilder.out._import(hfunc, ROOT.RooFit.RecycleConflictNodes())
        return self.modelBuilder.out.function(name)

    def doHistFuncFromXsecTools(self, higgs, quantity, varlist, production=None):
        # Translator mssm_xs_tools -> TH1D -> RooDataHist
        name  = self.quantity_map[quantity]['name']
        accesskey = self.quantity_map[quantity]['access']
        method = self.quantity_map[quantity]['method']
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
                value = getattr(self.mssm_inputs, method)(accesskey, x, y)
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
        # Computing scaling function for qqh contribution (SM-like Higgs) in context of MSSM
        # Assuming qqphi is already scaled to an appropriate SM 125.4 cross-section
        name  = "sf_qqphi_MSSM"

        accesskey_br = self.quantity_map['br']['access'].format(HIGGS=self.smlike)
        accesskey_br_SM = self.quantity_map['br_SM']['access']
        accesskey_vbf = self.quantity_map['xsec']['access'].format(PROD="vbf",HIGGS=self.smlike)
        accesskey_vbf_SM = self.quantity_map['xsec_SM']['access'].format(PROD="vbf", ADD="")
        accesskey_Wh = self.quantity_map['xsec']['access'].format(PROD="hs",HIGGS="W"+self.smlike)
        accesskey_Wh_SM = self.quantity_map['xsec_SM']['access'].format(PROD="hs",ADD="W")
        accesskey_Zh = self.quantity_map['xsec']['access'].format(PROD="hs",HIGGS="Z"+self.smlike)
        accesskey_Zh_SM = self.quantity_map['xsec_SM']['access'].format(PROD="hs",ADD="Z")

        print "Computing 'qqphi' scaling function from xsec tools"

        x_parname = varlist[0].GetName()
        x_binning = self.binning[self.scenario][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self.binning[self.scenario][y_parname]

        hist = ROOT.TH2D(name, name, len(x_binning)-1, x_binning, len(y_binning)-1, y_binning)
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):

                br_htautau = getattr(self.mssm_inputs, self.quantity_map['br']['method'])(accesskey_br, x, y)
                br_htautau_SM = getattr(self.mssm_inputs, self.quantity_map['br_SM']['method'])(accesskey_br_SM, x, y)
                xsec_vbf = getattr(self.mssm_inputs, self.quantity_map['xsec']['method'])(accesskey_vbf, x, y)
                xsec_vbf_SM = getattr(self.mssm_inputs, self.quantity_map['xsec_SM']['method'])(accesskey_vbf_SM, x, y)
                xsec_Wh = getattr(self.mssm_inputs, self.quantity_map['xsec']['method'])(accesskey_Wh, x, y)
                xsec_Wh_SM = getattr(self.mssm_inputs, self.quantity_map['xsec_SM']['method'])(accesskey_Wh_SM, x, y)
                xsec_Zh = getattr(self.mssm_inputs, self.quantity_map['xsec']['method'])(accesskey_Zh, x, y)
                xsec_Zh_SM = getattr(self.mssm_inputs, self.quantity_map['xsec_SM']['method'])(accesskey_Zh_SM, x, y)

                xsec = xsec_vbf + xsec_Wh + xsec_Zh
                xsec_SM = xsec_vbf_SM + xsec_Wh_SM + xsec_Zh_SM

                if br_htautau <= 0 and br_htautau_SM <= 0:
                    print "[WARNING]: Both BSM and SM BR predictions are <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting both to 1.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    br_htautau_SM = 1.
                    br_htautau = 1.
                elif br_htautau_SM <= 0:
                    print "[WARNING]: SM BR prediction is <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting to BSM prediction.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    br_htautau_SM = br_htautau

                if xsec <= 0 and xsec_SM <= 0:
                    print "[WARNING]: Both BSM and SM xsec predictions are <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting both to 1.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    xsec_SM = 1.
                    xsec = 1.
                elif br_htautau_SM <= 0:
                    print "[WARNING]: SM xsec prediction is <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting to BSM prediction.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    xsec_SM = xsec

                value = xsec / xsec_SM # xsec(mh) / xsec_SM(mh), correcting for mass dependence mh vs. 125.4 GeV
                value *= br_htautau / br_htautau_SM # br_htautau(mh) / br_htautau_SM(mh), correcting for mass dependence mh vs. 125.4 GeV
                value *= self.scaleforh # additional manual rescaling of light scalar h (default is 1.0)
                if self.use_hSM_difference:
                    value -= 1.0
                hist.SetBinContent(i_x+1, i_y+1, value)

        return self.doHistFunc(name, hist, varlist)

    def doHistFuncForGGH(self, varlist):
        # Computing scaling function for ggphi contribution (SM-like Higgs) in context of MSSM
        # Assuming ggphi is already scaled to an appropriate SM 125.4 cross-section
        name  = "sf_ggphi_MSSM"

        accesskey_xs = self.quantity_map['xsec']['access'].format(HIGGS=self.smlike,PROD='gg')
        accesskey_xs_SM = self.quantity_map['xsec_SM']['access'].format(ADD="",PROD='gg')
        accesskey_br = self.quantity_map['br']['access'].format(HIGGS=self.smlike)
        accesskey_br_SM = self.quantity_map['br_SM']['access']

        print "Computing 'ggphi' scaling function from xsec tools"

        x_parname = varlist[0].GetName()
        x_binning = self.binning[self.scenario][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self.binning[self.scenario][y_parname]

        hist = ROOT.TH2D(name, name, len(x_binning)-1, x_binning, len(y_binning)-1, y_binning)
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):

                xs_ggh = getattr(self.mssm_inputs, self.quantity_map['xsec']['method'])(accesskey_xs, x, y)
                xs_ggh_SM = getattr(self.mssm_inputs, self.quantity_map['xsec_SM']['method'])(accesskey_xs_SM, x, y)

                br_htautau = getattr(self.mssm_inputs, self.quantity_map['br']['method'])(accesskey_br, x, y)
                br_htautau_SM = getattr(self.mssm_inputs, self.quantity_map['br_SM']['method'])(accesskey_br_SM, x, y)

                if xs_ggh <= 0 and xs_ggh_SM <= 0:
                    print "[WARNING]: Both BSM and SM ggh xs predictions are <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting both to 1.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    xs_ggh_SM = 1.
                    xs_ggh = 1.
                elif xs_ggh_SM <= 0:
                    print "[WARNING]: SM ggh xs prediction is <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting to BSM prediction.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    xs_ggh_SM = xs_ggh

                if br_htautau <= 0 and br_htautau_SM <= 0:
                    print "[WARNING]: Both BSM and SM BR predictions are <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting both to 1.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    br_htautau_SM = 1.
                    br_htautau = 1.
                elif br_htautau_SM <= 0:
                    print "[WARNING]: SM BR prediction is <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting to BSM prediction.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    br_htautau_SM = br_htautau

                value =  xs_ggh / xs_ggh_SM * br_htautau / br_htautau_SM # xs(mh) * BR(mh) / (xs_SM(mh) * BR_SM(mh)) correcting for mass dependence mh vs. 125.4 GeV
                value *= self.scaleforh # additional manual rescaling of light scalar h (default is 1.0)
                hist.SetBinContent(i_x+1, i_y+1, value)

        return self.doHistFunc(name, hist, varlist)

    def doHistFuncForBBH(self, varlist):
        # Computing scaling function for bbphi contribution (SM-like Higgs) in context of MSSM
        # Assuming bbphi is **NOT** scaled to an appropriate SM 125.4 cross-section & BR, but to 1 pb
        name  = "sf_bbphi_MSSM"

        accesskey_xs = self.quantity_map['xsec']['access'].format(HIGGS=self.smlike,PROD='bb')
        accesskey_xs_SM = self.quantity_map['xsec_SM']['access'].format(ADD="",PROD='bb')
        accesskey_br = self.quantity_map['br']['access'].format(HIGGS=self.smlike)
        accesskey_br_SM = self.quantity_map['br_SM']['access']

        xs_bbh_SM125 = self.sm_predictions["xs_bb_SMH125"]
        br_htautau_SM125 = self.sm_predictions["br_SMH125_tautau"]

        print "Computing 'bbphi' scaling function from xsec tools"

        x_parname = varlist[0].GetName()
        x_binning = self.binning[self.scenario][x_parname]

        y_parname = varlist[1].GetName()
        y_binning = self.binning[self.scenario][y_parname]

        hist = ROOT.TH2D(name, name, len(x_binning)-1, x_binning, len(y_binning)-1, y_binning)
        for i_x, x in enumerate(x_binning):
            for i_y, y in enumerate(y_binning):

                xs_bbh = getattr(self.mssm_inputs, self.quantity_map['xsec']['method'])(accesskey_xs, x, y)
                xs_bbh_SM = getattr(self.mssm_inputs, self.quantity_map['xsec_SM']['method'])(accesskey_xs_SM, x, y)

                br_htautau = getattr(self.mssm_inputs, self.quantity_map['br']['method'])(accesskey_br, x, y)
                br_htautau_SM = getattr(self.mssm_inputs, self.quantity_map['br_SM']['method'])(accesskey_br_SM, x, y)

                if br_htautau <= 0 and br_htautau_SM <= 0:
                    print "[WARNING]: Both BSM and SM BR predictions are <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting both to 1.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    br_htautau_SM = 1.
                    br_htautau = 1.
                elif br_htautau_SM <= 0:
                    print "[WARNING]: SM BR prediction is <= 0 for {MASS}={MASSVAL}, tanb={TANBVAL}. Setting to BSM prediction.".format(MASS=self.massparameter, MASSVAL=x, TANBVAL=y)
                    br_htautau_SM = br_htautau

                # xs(mh) * (xs_SM(125.4)/xs_SM(mh)) * BR(mh) * (BR_SM(125.4)/BR_SM(mh)) correcting for mass dependence mh vs. 125.4 GeV
                value =  xs_bbh * (xs_bbh_SM125 / xs_bbh_SM) * br_htautau * (br_htautau_SM125 / br_htautau_SM)
                value *= self.scaleforh # additional manual rescaling of light scalar h (default is 1.0)
                hist.SetBinContent(i_x+1, i_y+1, value)

        return self.doHistFunc(name, hist, varlist)

    def doAsymPowSystematic(self, higgs, quantity, varlist, production, uncertainty):
        # Translator mssm_xs_tools -> TH1D -> RooDataHist -> Systematic
        name  = self.quantity_map[quantity]['name'].format(HIGGS=higgs, PROD=production)
        accesskey = self.quantity_map[quantity]['access'].format(HIGGS=higgs, PROD=production)
        uncertaintykey = self.uncertainty_map[production+uncertainty]
        method = self.quantity_map[quantity]['method']

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
                nominal  = getattr(self.mssm_inputs, method)(accesskey, x, y)
                value_hi = getattr(self.mssm_inputs, method)(accesskey+uncertaintykey.format(VAR='up'), x, y)
                value_lo = getattr(self.mssm_inputs, method)(accesskey+uncertaintykey.format(VAR='down'), x, y)
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

        fractions_sm = ROOT.TFile.Open(self.ggHatNLO, 'read')
        w_sm = fractions_sm.Get("w")

        template_X = X
        if self.scenario == "mh1125_CPV":
            template_X = self.cpv_template_pairs[X]
            w_sm.var("m{HIGGS}".format(HIGGS=self.cpv_template_pairs[X])).SetName("m{HIGGS}".format(HIGGS=X))

        for loopcontrib in ['t','b','i']:
            func = w_sm.function('gg{X}_{LC}_MSSM_frac'.format(X=template_X, LC=loopcontrib))
            getattr(self.modelBuilder.out, 'import')(func, ROOT.RooFit.RecycleConflictNodes())
            self.modelBuilder.out.factory('prod::%s(%s,%s)' % (name.format(X=X, LC="_"+loopcontrib), name.format(X=X, LC=""), "gg%s_%s_MSSM_frac" % (template_X,loopcontrib))) #multiply t,b,i fractions with xsec at NNLO


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
        self.modelBuilder.out.var(self.massparameter).setConstant(True) # either mA or mHp
        self.modelBuilder.out.var('tanb').setConstant(True)

        bsm_proc_match = "(gg(A|H|h|H3|H2|H1)_(t|i|b)|bb(A|H|h|H3|H2|H1))(|_lowmass)"
        if self.replace_with_sm125:
            bsm_proc_match = "(gg(A|H|h|H3|H2|H1)_(t|i|b)|bb(A|{BSMSCALAR}|H3|H2))(|_lowmass)".format(BSMSCALAR=self.bsmscalar).replace("||","|") # need the last fix in case BSMSCALAR=""

        for proc in self.PROC_SETS:
            terms = []
            X = proc.split('_')[0].replace('gg','').replace('bb','')
            ggphi_smlike_match = re.match('gg{SMLIKE}$'.format(SMLIKE=self.smlike), proc)
            bbphi_smlike_match = re.match('bb{SMLIKE}$'.format(SMLIKE=self.smlike), proc)
            if "H125" in proc: # cover SM H125 processes first
                terms = [self.sigNorms[False]]
            elif re.match(bsm_proc_match, proc): # not SM-like BSMSCALAR: either h or H
                terms = ['xs_%s' %proc.replace("_lowmass",""), 'br_%stautau'%X]
                terms += ['r']
                terms += [self.sigNorms[True]]
            elif re.match('(qq{SMLIKE}|Z{SMLIKE}|W{SMLIKE})$'.format(SMLIKE=self.smlike), proc): # always done
                terms = [self.sigNorms[True], 'r', 'sf_qqphi_MSSM']
            elif ggphi_smlike_match: # always done
                terms = [self.sigNorms[True], 'r', 'sf_ggphi_MSSM']
            elif bbphi_smlike_match: # considered, in case it is not in the second 'if' case with bsm_proc_match
                terms = [self.sigNorms[True], 'r', 'sf_bbphi_MSSM']

            if self.scenario == "mh1125_CPV" and X in ['H2', 'H3']:
                for xx in ['bb', 'gg']:
                    if xx in proc:
                        terms.append('expr::interference_{PROD}_{HIGGS}(\"1.0 + @0\", int_{PROD}{HIGGS}_tautau)'.format(PROD=xx, HIGGS=X))

            # Now scan terms and add theory uncerts
            extra = []
            for term in terms:
                if term in self.SYST_DICT:
                    extra += self.SYST_DICT[term]
            terms += extra
            if ggphi_smlike_match and self.use_hSM_difference:
                self.modelBuilder.factory_('prod::bsm_scaling_%s(%s)'%(proc,','.join(terms))) # Add scaling of BSM process: mu*SF = x*r*SF
                self.modelBuilder.factory_('expr::scaling_%s(\"(@0 - @1 * @2)\", %s)'%(proc,','.join(["bsm_scaling_%s"%proc,"x","r"]))) # Add difference between BSM prediction and mu*SM prediction.
            elif bbphi_smlike_match and self.use_hSM_difference and self.replace_with_sm125:
                self.modelBuilder.factory_('prod::bsm_scaling_%s(%s)'%(proc,','.join(terms)))
                self.modelBuilder.factory_('expr::scaling_%s(\"(@0 - @1 * @2 * %s * %s)\", %s)'%(proc,str(self.sm_predictions["xs_bb_SMH125"]),str(self.sm_predictions["br_SMH125_tautau"]),','.join(["bsm_scaling_%s"%proc,"x","r"])))
            else:
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
        mass = ROOT.RooRealVar(self.massparameter, 'm_{A} [GeV]' if self.massparameter == 'mA' else 'm_{H^{+}} [GeV]', 160.) # the other case would be 'mHp'
        tanb = ROOT.RooRealVar('tanb', 'tan#beta', 5.5)
        pars = [mass, tanb]

        self.mssm_inputs = mssm_xs_tools(self.filename, False, 1) # syntax: model filename, Flag for interpolation ('True' or 'False'), verbosity level

        # qqphi, Zphi and Wphi  added always in this setup
        self.doHistFuncForQQH(pars)
        self.PROC_SETS.extend(['qq'+self.smlike, 'Z'+self.smlike, 'W'+self.smlike])

        # adding ggphi & bbphi as 125 templates only if requested
        if self.replace_with_sm125:

            self.doHistFuncForGGH(pars)
            self.PROC_SETS.append('gg'+self.smlike)

            self.doHistFuncForBBH(pars) # no need to add since added later

        procs = ['H1', 'H2', 'H3'] if self.scenario == "mh1125_CPV" else ['h', 'H', 'A']

        for X in procs:
            if self.massparameter.replace('m','') == X: # don't create histogram for 'A' in cases, where its mass is a model-parameter
                continue
            self.doHistFuncFromXsecTools(X, "mass", pars) # syntax: Higgs-Boson, mass attribute, parameters

        for X in procs:
            if self.scenario == "mh1125_CPV":
                self.doHistFuncFromXsecTools(X, "interference", pars, production="gg")
                self.doHistFuncFromXsecTools(X, "interference", pars, production="bb")
            else:
                self.doHistFuncFromXsecTools(X, "yukawa_top", pars)
                self.doHistFuncFromXsecTools(X, "yukawa_bottom", pars)

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

            for loopcontrib in ['t','b','i']:
                self.SYST_DICT['xs_gg%s_%s' % (X, loopcontrib)].append('systeff_xs_gg%s_MSSM_scale' %X)
                self.SYST_DICT['xs_gg%s_%s' % (X, loopcontrib)].append('systeff_xs_gg%s_MSSM_pdfas' %X)

            self.SYST_DICT['xs_bb%s' %X].append('systeff_xs_bb%s_MSSM_total' %X)

            # Make a note of what we've built, will be used to create scaling expressions later
            self.PROC_SETS.append('bb%s'%X)
            if X != self.smlike:
                self.PROC_SETS.append('bb%s_lowmass'%X)
            self.PROC_SETS.extend(['gg%s_t'%X, 'gg%s_b'%X, 'gg%s_i'%X, 'gg%s_t_lowmass'%X, 'gg%s_b_lowmass'%X, 'gg%s_i_lowmass'%X])

        # Add BSM systematic also in case SM125 templates are used for ggphi and bbphi
        if self.replace_with_sm125:
             self.SYST_DICT["sf_ggphi_MSSM"].append('systeff_xs_gg%s_MSSM_scale' %self.smlike)
             self.SYST_DICT["sf_ggphi_MSSM"].append('systeff_xs_gg%s_MSSM_pdfas' %self.smlike)
             self.SYST_DICT["sf_bbphi_MSSM"].append('systeff_xs_bb%s_MSSM_total' %self.smlike)

        # And the SM terms
        self.PROC_SETS.extend(['ggH125', 'qqH125', 'ZH125', 'WH125', 'bbH125'])
        if self.debug_output:
            self.debug_output.Close()


MSSMvsSM = MSSMvsSMHiggsModel()
