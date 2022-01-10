import ROOT
import math
import argparse
import os
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

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname,extraFFuncerts=False,year=2018):
    directory = infile.Get(dirname)
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
 
        if not isinstance(histo,ROOT.TDirectory):
            WriteToTFile(histo, outfile, write_dirname+"/"+key.GetName())
            if histo.GetName() == 'jetFakes' and extraFFuncerts:
              nbins = histo.GetNbinsX()
              x=4 # number of pt_tt bins
              nxbins = nbins/x
              for i in range(1,5):
                hup = histo.Clone()
                hdown = histo.Clone()
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
newfilename=filename.replace('.root','_origional.root')

os.system('mv %(filename)s %(newfilename)s'%vars())

filename_2 = filename.replace('mt_tot_puppi','m_sv_puppi')
filename_3 = filename.replace('mt_tot_puppi','m_sv_VS_pt_tt')

original_file = ROOT.TFile(newfilename)
output_file_1 = ROOT.TFile(filename,"RECREATE")
output_file_2 = ROOT.TFile(filename_2,"RECREATE")
output_file_3 = ROOT.TFile(filename_3,"RECREATE")

match_2 = '_svfit'
match_3 = '_pt_tt_vs_svfit'

year=''
if '2018' in filename_3.split('/')[-1]: year = 2018
if '2017' in filename_3.split('/')[-1]: year = 2017
if '2016' in filename_3.split('/')[-1]: year = 2016

for key in original_file.GetListOfKeys():
    print key
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        print dirname
        #directory = original_file.Get(dirname)

        new_out_name=None
        new_dir_name=None

        if match_2 not in dirname and match_3 not in dirname:
          new_out_name = output_file_1
          new_dir_name = dirname

        if match_2 in dirname and match_3 not in dirname:
          new_out_name = output_file_2
          new_dir_name = dirname.replace(match_2,'')

        if match_3 in dirname:
          new_out_name = output_file_3
          new_dir_name = dirname.replace(match_3,'')
 
        if new_out_name: getHistogramAndWriteToFile(original_file,new_out_name,dirname,new_dir_name,match_3,year)
 
