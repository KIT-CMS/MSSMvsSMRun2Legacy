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
        if isinstance(histo,ROOT.TH1F) or isinstance(histo,ROOT.TH1D):
            histonew=ROOT.TH1D(histo.GetName(), histo.GetTitle(),histo.GetNbinsX(),0,histo.GetNbinsX())
            for i in range(0,histo.GetNbinsX()+2):
              c = histo.GetBinContent(i)
              e = histo.GetBinError(i)
              histonew.SetBinContent(i,c)
              histonew.SetBinError(i,e)
            WriteToTFile(histonew, outfile, write_dirname+"/"+key.GetName())

parser = argparse.ArgumentParser()
parser.add_argument('--file', '-f', help= 'File from which subdirectories need to be dropped')
args = parser.parse_args()
filename = args.file
newfilename=filename.replace('.root','_unmodified.root')

os.system('mv %(filename)s %(newfilename)s'%vars())

original_file = ROOT.TFile(newfilename)
output_file = ROOT.TFile(filename,"RECREATE")

for key in original_file.GetListOfKeys():
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        convertHistogramAndWriteToFile(original_file,output_file,key.GetName(),dirname)
