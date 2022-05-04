import ROOT
import array

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bOnly', help= 'Use b-only fit result', action='store_true')
args = parser.parse_args()


if args.bOnly: fout = ROOT.TFile('shapes_cbyears_bOnly.root','RECREATE')
else: fout = ROOT.TFile('shapes_cbyears.root','RECREATE')

cb_procs = ['TotalBkg',	'TotalProcs', 'TotalSig', 'data_obs']

procs = {}
procs['lt'] = ['EMB','TTL','VVL','ZL','bbH125','bbh','ggH125','ggh_b','ggh_i','ggh_t','jetFakes','qqH125']
procs['em'] = ['EMB','QCD','TTL','VVL','W','WHWW125','ZHWW125','ZL','bbH125','bbh','ggH125','ggHWW125','ggh_b','ggh_i','ggh_t','qqH125','qqHWW125']
procs['tt'] = ['EMB','TTL','VVL','ZL','bbH125','bbh','ggH125','ggh_b','ggh_i','ggh_t','jetFakes','qqH125','wFakes']

dir_map = {
  'tt_32' : 'tt_Nbtag0',
  'tt_35' : 'tt_NbtagGt1',
  'mt_32' : 'mt_Nbtag0_MTLt40',
  'mt_35' : 'mt_NbtagGt1_MTLt40',
  'et_32' : 'et_Nbtag0_MTLt40',
  'et_35' : 'et_NbtagGt1_MTLt40',
  'em_33' : 'em_Nbtag0_DZetam10To30',
  'em_32' : 'em_Nbtag0_DZetaGt30',
  'em_36' : 'em_NbtagGt1_DZetam10To30',
  'em_35' : 'em_NbtagGt1_DZetaGt30',
}


for c in ['lt','tt','em']:
  #bins = [32,35,432,332]
  bins = [35,432,332,232,132]
  if args.bOnly: 
    bins+=['32_mt_tot','35_mt_tot']
    if c in ['lt','em']: bins+=['33_mt_tot','36_mt_tot']
    if c in ['em']: bins+=['34_mt_tot','37_mt_tot']
  for b in bins:

    out_dir = 'htt_%(c)s_%(b)s_postfit' % vars()
    fout.mkdir(out_dir)

    if args.bOnly: fin = ROOT.TFile('shapes_cbyears_bOnly_%(c)s_%(b)s.root' % vars())
    else: fin = ROOT.TFile('shapes_cbyears_%(c)s_%(b)s.root' % vars())

    print fin, 'shapes_cbyears_bOnly_%(c)s_%(b)s.root' % vars()


    if c == 'lt': chans = ['mt','et']
    else: chans = [c]

    # get totals first
    for x in cb_procs:
      h = fin.Get('postfit/%(x)s' % vars())
      fout.cd(out_dir)
      h.Write()


    # loop over each process, hadd histograms for different years / channels then write the total to the file
    for x in procs[c]:
      h = fin.Get('postfit/data_obs').Clone()
      h.Reset()
      h.SetName(x)
      for y in [2016, 2017, 2018]:
        for chan in chans:
           indir = 'htt_%(chan)s_%(b)s_%(y)s_postfit' % vars()
           htemp = fin.Get('%(indir)s/%(x)s' % vars())
           if chan =='em':
             if isinstance(b,str) and 'mt_tot' in b: b_= '%i_mt_tot' % (int(b.split('_')[0])+1)
             else: b_=b+1
             indir2 = 'htt_%(chan)s_%(b_)s_%(y)s_postfit' % vars()
             htemp2 = fin.Get('%(indir2)s/%(x)s' % vars())
             if isinstance(htemp2,ROOT.TH1D) or isinstance(htemp2,ROOT.TH1F): htemp.Add(htemp2)
           if isinstance(htemp,ROOT.TH1D) or isinstance(htemp,ROOT.TH1F): h.Add(htemp)

      fout.cd(out_dir)
      h.Write()

    # for high mass categories we also add 1 TeV VLQ signal scales to bestfit point (gU=1.2)
    if b not in [33,36] or True: continue
    h1 = fin.Get('postfit/data_obs').Clone()
    h2 = fin.Get('postfit/data_obs').Clone()
    h1.Reset()
    h2.Reset()
    h1.SetName('VLQ_s')
    h2.SetName('VLQ_i')
    new_bins=[]
    for i in range(1,h1.GetNbinsX()+2):
      new_bins.append(h1.GetBinLowEdge(i))
    new_bins= array.array('d',new_bins)
    for y in [2016, 2017, 2018]:
      for chan in chans:
         print 'shapes/%(y)s/%(chan)s/vlq.inputs-mssm-vs-sm-Run%(y)s-mt_tot_puppi.root' % vars()
         fin_vlq = ROOT.TFile('shapes/%(y)s/%(chan)s/vlq.inputs-mssm-vs-sm-Run%(y)s-mt_tot_puppi.root' % vars())
         indir = dir_map['%(chan)s_%(b)s' % vars()]
         htemp1 = fin_vlq.Get('%(indir)s/VLQ_betaRd33_0_matched_M_1000' % vars())
         htemp2 = fin_vlq.Get('%(indir)s/VLQ_betaRd33_0_matched_interference_M_1000' % vars())
         print '%(indir)s/VLQ_betaRd33_0_matched_M_1000' % vars()
         print htemp1
         htemp1 = htemp1.Rebin(len(new_bins)-1,"",new_bins)
         htemp2 = htemp2.Rebin(len(new_bins)-1,"",new_bins)
         h1.Add(htemp1)
         h2.Add(htemp2)

    fout.cd(out_dir)

    # scale vlq to gU=1.2
    gU=1.2
    h1.Scale(1.2**4)
    h2.Scale(1.2**2)

    # remove stat fluctuations in low mt_tot bins
    for i in range(1,h1.FindBin(100)):
      c1 = h1.GetBinContent(i)   
      c2 = h1.GetBinContent(i)
      if c1+c2>0: 
        h1.SetBinContent(i,0.)   
        h2.SetBinContent(i,0.)   

    h1.Write()
    h2.Write()

          

    

