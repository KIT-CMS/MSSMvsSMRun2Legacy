#!/usr/bin/env python
#This script drops unnecessary subdirectories from the shape files - the original
#input file is saved under the same name but with _full_file appended
import ROOT
import math
import argparse
import json
import sys
import os
import fnmatch
ROOT.TH1.AddDirectory(False)

def WriteToTFile(obj, file, path):
    file.cd()
    as_vec = path.split('/')
    if len(as_vec) >= 1:
        for i in xrange(0, len(as_vec)-1):
            if not ROOT.gDirectory.GetDirectory(as_vec[i]):
                ROOT.gDirectory.mkdir(as_vec[i])
            ROOT.gDirectory.cd(as_vec[i])
    if not ROOT.gDirectory.FindKey(as_vec[-1]):
        obj.SetName(as_vec[-1])
        ROOT.gDirectory.WriteTObject(obj, as_vec[-1])
    ROOT.gDirectory.cd('/')

def ScaleJetFakes(histo,chan, mt_cat, year):

  func_map = {
               'mt': {2016: {'mTloose': '(0.987985+-0.00232843*x+2.10697e-05*pow(x,2)+-4.30698e-08*pow(x,3))', 'mTtight': '(0.919909+-0.00175173*x+2.10265e-05*pow(x,2)+-4.33267e-08*pow(x,3))'}, 2017: {'mTloose': '(0.890062+-0.00208239*x+1.95682e-05*pow(x,2)+-4.13481e-08*pow(x,3))', 'mTtight': '(0.932662+-0.00324539*x+2.64541e-05*pow(x,2)+-4.83053e-08*pow(x,3))'}, 2018: {'mTloose': '(0.881833+-0.00217206*x+1.97969e-05*pow(x,2)+-3.80306e-08*pow(x,3))', 'mTtight': '(0.862101+-0.000914674*x+9.13855e-06*pow(x,2)+-1.49594e-08*pow(x,3))'}},
                'et': {2016: {'mTloose': '(0.844816+-0.000543283*x+9.73709e-06*pow(x,2)+-1.2664e-08*pow(x,3))', 'mTtight': '(0.70624+0.00191175*x+-2.66422e-06*pow(x,2)+5.20671e-09*pow(x,3))'}, 2017: {'mTloose': '(0.930003+-0.0026238*x+2.06055e-05*pow(x,2)+-3.554e-08*pow(x,3))', 'mTtight': '(1.05848+-0.00450755*x+2.98441e-05*pow(x,2)+-5.30103e-08*pow(x,3))'}, 2018: {'mTloose': '(0.938801+-0.0032487*x+2.69341e-05*pow(x,2)+-5.57956e-08*pow(x,3))', 'mTtight': '(0.748626+0.00146163*x+-6.11159e-06*pow(x,2)+1.222e-08*pow(x,3))'}}
             }

  print histo.GetName(), dirname, chan, mt_cat, year, func_map[chan][year][mt_cat]
  func = ROOT.TF1('func',func_map[chan][year][mt_cat])  

  histonew=histo.Clone()
  for i in range(1,histonew.GetNbinsX()+1):
    b = histo.GetBinCenter(i)
    sf = func.Eval(b)
    c = histo.GetBinContent(i)
    e = histo.GetBinError(i)
    histonew.SetBinContent(i,c*sf)
    histonew.SetBinError(i,e*sf)
    
  return histonew

def convertHistogramAndWriteToFile(infile,outfile,dirname,write_dirname,extraFFuncerts=False,year=2018):
    directory = infile.Get(dirname)
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1F) or isinstance(histo,ROOT.TH1D):
            if 'jetFakes' in histo.GetName() and ('NbtagGt1' in dirname):
             if dirname.startswith('mt'): chan='mt' 
             if dirname.startswith('et'): chan='et'
             if 'MTLt40' in dirname: mt_cat = 'mTtight'
             if 'MT40To70' in dirname: mt_cat = 'mTloose' 
             histonew = ScaleJetFakes(histo,chan, mt_cat, year)
             WriteToTFile(histonew, outfile, write_dirname+"/"+key.GetName())
             if histo.GetName() == 'jetFakes': 
               new_name='jetFakes_CMS_ff_total_ttbar_msv_shape_syst_%(mt_cat)s_%(year)s' % vars()
               WriteToTFile(histo, outfile, write_dirname+"/"+new_name+'Down') 
               histonewup = ScaleJetFakes(histonew,chan, mt_cat, year)
               WriteToTFile(histonewup, outfile, write_dirname+"/"+new_name+'Up') 
            else: WriteToTFile(histo, outfile, write_dirname+"/"+key.GetName())

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which subdirectories need to be dropped')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','_unmodified.root')

os.system('mv %(filename)s %(newfilename)s'%vars())

original_file = ROOT.TFile(newfilename)
output_file = ROOT.TFile(filename,"RECREATE")

year=''
if '2018' in filename.split('/')[-1]: year = 2018
if '2017' in filename.split('/')[-1]: year = 2017
if '2016' in filename.split('/')[-1]: year = 2016

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        convertHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname,True,year)
