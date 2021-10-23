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

def getHistogramAndWriteToFile(infile,outfile,dirname,write_dirname):
    directory = infile.Get(dirname)
    for key in directory.GetListOfKeys():
        histo = directory.Get(key.GetName())
 
        if not isinstance(histo,ROOT.TDirectory):
            WriteToTFile(histo, outfile, write_dirname+"/"+key.GetName())

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

for key in original_file.GetListOfKeys():
    print key
    if isinstance(original_file.Get(key.GetName()),ROOT.TDirectory):
        dirname=key.GetName()
        print dirname
        #directory = original_file.Get(dirname)

        new_out_name=None
        new_out_name=None

        if match_2 not in dirname and match_3 not in dirname:
          new_out_name = output_file_1
          new_dir_name = dirname
          print 'test1' 

        if match_2 in dirname and match_3 not in dirname:
          new_out_name = output_file_2
          new_dir_name = dirname.replace(match_2,'')
          print 'test2' 

        if match_3 in dirname:
          new_out_name = output_file_3
          new_dir_name = dirname.replace(match_3,'')
          print 'test3' 
 
        if new_out_name: getHistogramAndWriteToFile(original_file,new_out_name,dirname,new_dir_name)
 
