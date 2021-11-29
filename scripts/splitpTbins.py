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

def convertHistogramAndWriteToFile(infile,outfile,dirname,write_dirname):
    directory = infile.Get(dirname)
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
        if histo.GetNbinsX()!=120: continue
        if isinstance(histo,ROOT.TH1F) or isinstance(histo,ROOT.TH1D):

            histo_template=ROOT.TH1D(histo.GetName(), histo.GetTitle(),30,0,300)
            for i in range (0,4):
              histonew = histo_template.Clone()
              for j in range(1,31):
		c = histo.GetBinContent(j+i*30)
                e = histo.GetBinError(j+i*30)
                histonew.SetBinContent(j,c)
                histonew.SetBinError(j,e)

              if i==0: write_dirname_extra = write_dirname+'_pT_0To50'
              if i==1: write_dirname_extra = write_dirname+'_pT_50To100'
              if i==2: write_dirname_extra = write_dirname+'_pT_100To200'
              if i==3: write_dirname_extra = write_dirname+'_pT_GT200'

              WriteToTFile(histonew, outfile, write_dirname_extra+"/"+key.GetName())

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which subdirectories need to be dropped')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','_splitpT.root')

original_file = ROOT.TFile(filename)
output_file = ROOT.TFile(newfilename,"RECREATE")

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        convertHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname)
