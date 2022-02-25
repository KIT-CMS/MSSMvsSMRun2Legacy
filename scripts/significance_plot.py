import ROOT
from array import array
import math
from collections import OrderedDict
import CombineHarvester.CombineTools.plotting as plot

plot.ModTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.06)

def DrawTitle(pad, text, align, scale=1):
    pad_backup = ROOT.gPad
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()

    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if pad_ratio < 1.:
        pad_ratio = 1.

    textSize = 0.6
    textOffset = 0.2

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(textSize * t * pad_ratio * scale)

    y_off = 1 - t + textOffset * t
    if align == 1:
        latex.SetTextAlign(11)
    if align == 1:
        latex.DrawLatex(l, y_off, text)
    if align == 2:
        latex.SetTextAlign(21)
    if align == 2:
        latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text)
    if align == 3:
        latex.SetTextAlign(31)
    if align == 3:
        latex.DrawLatex(1 - r, y_off, text)
    pad_backup.cd()

scenario = OrderedDict()

scenario["mt_tot"] = "m_{T}^{tot}"
scenario["lowmass"] = "m_{#tau#tau} vs p_{T}^{#tau#tau}"

colour = {
           "mt_tot":4,
           "lowmass":4,
           }

for proc in ['gg','bb']:

  c = ROOT.TCanvas('c','c',700,700)
  plot.Set(c, Tickx=1, Ticky=1)
  c.SetLogy()
  c.SetLogx()
  
  limit_dict = {}
  first_pass = True
  for key, val in scenario.items():
     
    # change so in bin center
    
    limit_dict[key] = {} 
  
    limit_dict[key]["p_value"] = ROOT.TGraph()
    limit_dict[key]["significance"] = ROOT.TGraph()
    
    
    i = 1
    if key == 'lowmass':
      masses = [60,80,95,100,120,125,130,140,160,180,200]
    else: 
      masses = [250,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2300,2600,2900,3200,3500]
    for m in masses:
      if key == 'lowmass':
        f = ROOT.TFile("model_independent_limits/Feb07_all_all_bsm-model-indep/combined/cmb/higgsCombine.%(proc)sH.v2.Significance.mH%(m)s.root" % vars())
      else: 
        f = ROOT.TFile("model_independent_limits/Jan12_mt_tot_all_all_bsm-model-indep/combined/cmb/higgsCombine.%(proc)sH.Significance.mH%(m)s.root" % vars())
      t = f.Get('limit')
      for event in t:
        limit_dict[key]["significance"].SetPoint(limit_dict[key]["significance"].GetN(),m,event.limit)
      for event in t:
        limit_dict[key]["p_value"].SetPoint(limit_dict[key]["p_value"].GetN(),m,ROOT.Math.normal_cdf_c(event.limit))
      i+=1
    
    #limit_dict[key]["p_value"].Print("all")
    print key
    limit_dict[key]["significance"].Print("all")
  
    if first_pass:
  
     
      limit_dict[key]["p_value"].GetXaxis().SetLimits(50,3600) 
      limit_dict[key]["p_value"].GetXaxis().SetRangeUser(50,3600)
      limit_dict[key]["p_value"].Draw("apl")
      limit_dict[key]["p_value"].SetLineWidth(2)
      #limit_dict[key]["p_value"].GetXaxis().SetTitleSize(0.04)
      limit_dict[key]["p_value"].GetXaxis().SetTitle('m_{#phi} (GeV)')
      limit_dict[key]["p_value"].GetYaxis().SetTitle('Local p-value')
#      limit_dict[key]["p_value"].GetYaxis().SetTitleOffset(1.2)
      limit_dict[key]["p_value"].GetXaxis().SetNoExponent()
      limit_dict[key]["p_value"].GetXaxis().SetMoreLogLabels()
      limit_dict[key]["p_value"].GetXaxis().SetLabelOffset(limit_dict[key]["p_value"].GetXaxis().GetLabelOffset()*2)
      limit_dict[key]["p_value"].SetLineColor(colour[key])
      limit_dict[key]["p_value"].SetMarkerColor(colour[key])
      limit_dict[key]["p_value"].SetMarkerStyle(15)
      limit_dict[key]["p_value"].SetMaximum(1)
      limit_dict[key]["p_value"].SetMinimum(0.0004)
  
      l = ROOT.TLegend(0.40,0.15,0.87,0.4)
      l.SetFillStyle(0)
      l.SetBorderSize(0)
      l.SetTextSize(0.04)
      l.AddEntry(limit_dict[key]["p_value"],val,'l')
  
      first_pass = False
    else:
      limit_dict[key]["p_value"].SetLineWidth(2)
      limit_dict[key]["p_value"].SetLineColor(colour[key])
      limit_dict[key]["p_value"].SetMarkerColor(colour[key])
      limit_dict[key]["p_value"].SetMarkerStyle(15)
      limit_dict[key]["p_value"].GetXaxis().SetRangeUser(50,3600)
      limit_dict[key]["p_value"].Draw("plsame")
  
      c.Update()
      l.AddEntry(limit_dict[key]["p_value"],val,'lp')
  
  
  #l.Draw()
  
  y = [ROOT.Math.normal_cdf_c(1.0),ROOT.Math.normal_cdf_c(2.0),ROOT.Math.normal_cdf_c(3.0)]
  line = []
  latex = []
  for i in range(0,len(y)):
    line.append(ROOT.TLine(50.,y[i],3600,y[i]))
    line[i].SetLineWidth(2)
    line[i].SetLineStyle(2)
    line[i].SetLineColor(13)
    line[i].Draw()
  
    latex.append(ROOT.TLatex())
    #latex[i].SetNDC()
    latex[i].SetTextAngle(0)
    latex[i].SetTextAlign(12)
    latex[i].SetTextSize(0.04)
    latex[i].SetTextColor(13)
    t = c.GetTopMargin()
    b = c.GetBottomMargin()
    r = c.GetRightMargin()
    #latex[i].DrawLatex(1.01-r,-b+((1-math.log(y[i])/math.log(0.001)*(1-t-b))),str(i+1) + "#sigma")
    latex[i].DrawLatex(3800.,y[i],str(i+1) + "#sigma")
  

  plot.DrawTitle(c, '138 fb^{-1} (13 TeV)', 3)
  plot.DrawCMSLogo(c, 'CMS', '', 1, 0.045, 0.05, 1.0, '', 0.9)
 
  line2 =  ROOT.TLine(225.,limit_dict["lowmass"]["p_value"].GetMinimum(),225,1.)
  line2.Draw()
  latex2 = ROOT.TLatex()
  latex2.SetNDC()
  latex2.SetTextAngle(0)
  latex2.SetTextAlign(12)
  latex2.SetTextFont(42)
  latex2.SetTextSize(0.04)
  latex2.DrawLatex(0.19,0.9, 'Low-mass')
  latex2.DrawLatex(0.45,0.9, 'High-mass')

  latex2.SetTextSize(0.05)
  latex2.DrawLatex(0.8,0.16, '%(proc)s#phi' % vars())
  c.Print('significance_plot_%(proc)sH.pdf' % vars())
