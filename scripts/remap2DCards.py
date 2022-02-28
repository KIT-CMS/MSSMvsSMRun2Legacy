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

def convertHistogramAndWriteToFile(infile,outfile,dirname,write_dirname,extraFFuncerts=False,year=2018):
    directory = infile.Get(dirname)
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if isinstance(histo,ROOT.TH1F) or isinstance(histo,ROOT.TH1D):
            histonew=ROOT.TH1D(histo.GetName(), histo.GetTitle(),histo.GetNbinsX(),0,histo.GetNbinsX())
            for i in range(0,histo.GetNbinsX()+2):
              c = histo.GetBinContent(i)
              e = histo.GetBinError(i)
              histonew.SetBinContent(i,c)
              histonew.SetBinError(i,e)
            WriteToTFile(histonew, outfile, write_dirname+"/"+key.GetName())

            if histonew.GetName() == 'jetFakes' and extraFFuncerts:
              nbins = histonew.GetNbinsX()
              x=4 # number of pt_tt bins
              nxbins = nbins/x
              for i in range(1,5):
                hup = histonew.Clone()
                hdown = histonew.Clone()
                for j in range(1+(i-1)*nxbins,1+i*nxbins):
                  hup.SetBinContent(j,hup.GetBinContent(j)*1.05)
                  hdown.SetBinContent(j,hdown.GetBinContent(j)*0.95)
                new_name='jetFakes_CMS_ff_total_syst_pt_tt_bin%(i)i_%(write_dirname)s_%(year)s' % vars()
                hup.SetName(new_name+'Up')
                hdown.SetName(new_name+'Down')
                WriteToTFile(hup, outfile, write_dirname+"/"+new_name+'Up')
                WriteToTFile(hdown, outfile, write_dirname+"/"+new_name+'Down')

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
        convertHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname,'m_sv_VS_pt_tt' in newfilename,year)
