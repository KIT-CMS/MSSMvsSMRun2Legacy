import CombineHarvester.CombineTools.plotting as plot
import ROOT
import re
import math
import argparse
import json
import numpy as np
import sys
import os
import fnmatch
from array import array

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(False)

#example commands:

#unrolled prefit plot for nobtag category with blinded bins
#python scripts/postFitPlotJetFakes.py --mode prefit --file_dir htt_mt_32_2018 -f shapes_mt_32_2018.root --ratio  --log_y --manual_blind
#1D prefit plot for btag category with blinded bins
#python scripts/postFitPlotJetFakes.py --mode prefit --file_dir htt_mt_35_2018 -f shapes_mt_35_2018.root --ratio  --manual_blind

def DrawCMSLogo(pad, cmsText, extraText, iPosX, relPosX, relPosY, relExtraDY, extraText2='', cmsTextSize=0.8,relExtraDX=0.0):
    """Blah

    Args:
        pad (TYPE): Description
        cmsText (TYPE): Description
        extraText (TYPE): Description
        iPosX (TYPE): Description
        relPosX (TYPE): Description
        relPosY (TYPE): Description
        relExtraDY (TYPE): Description
        extraText2 (str): Description
        cmsTextSize (float): Description

    Returns:
        TYPE: Description
    """
    pad.cd()
    cmsTextFont = 62  # default is helvetic-bold

    writeExtraText = len(extraText) > 0
    writeExtraText2 = len(extraText2) > 0
    extraTextFont = 52

    # text sizes and text offsets with respect to the top frame
    # in unit of the top margin size
    lumiTextOffset = 0.2
    # cmsTextSize = 0.8
    # float cmsTextOffset    = 0.1;  // only used in outOfFrame version

    # ratio of 'CMS' and extra text size
    extraOverCmsTextSize = 0.76

    outOfFrame = False
    if iPosX / 10 == 0:
        outOfFrame = True

    alignY_ = 3
    alignX_ = 2
    if (iPosX / 10 == 0):
        alignX_ = 1
    if (iPosX == 0):
        alignX_ = 1
    if (iPosX == 0):
        alignY_ = 1
    if (iPosX / 10 == 1):
        alignX_ = 1
    if (iPosX / 10 == 2):
        alignX_ = 2
    if (iPosX / 10 == 3):
        alignX_ = 3
    # if (iPosX == 0): relPosX = 0.14
    align_ = 10 * alignX_ + alignY_

    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    extraTextSize = extraOverCmsTextSize * cmsTextSize
    pad_ratio = (float(pad.GetWh()) * pad.GetAbsHNDC()) / \
        (float(pad.GetWw()) * pad.GetAbsWNDC())
    if (pad_ratio < 1.):
        pad_ratio = 1.

    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText)

    posX_ = 0
    if iPosX % 10 <= 1:
        posX_ = l + relPosX * (1 - l - r)
    elif (iPosX % 10 == 2):
        posX_ = l + 0.5 * (1 - l - r)
    elif (iPosX % 10 == 3):
        posX_ = 1 - r - relPosX * (1 - l - r)

    posY_ = 1 - t - relPosY * (1 - t - b)
    if not outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextSize(cmsTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, cmsText)
        if writeExtraText:
            latex.SetTextFont(extraTextFont)
            latex.SetTextAlign(align_)
            latex.SetTextSize(extraTextSize * t * pad_ratio)
            latex.DrawLatex(
                posX_- relExtraDX, posY_ - relExtraDY * cmsTextSize * t, extraText)
            if writeExtraText2:
                latex.DrawLatex(
                    posX_, posY_ - 1.8 * relExtraDY * cmsTextSize * t, extraText2)
    elif writeExtraText:
        if iPosX == 0:
            posX_ = l + relPosX * (1 - l - r) + relExtraDX
            posY_ = 1 - t + lumiTextOffset * t
        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(extraTextSize * t * pad_ratio)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_, posY_, extraText)


def getHistogram(fname, histname, dirname='', postfitmode='prefit', allowEmpty=False, logx=False):
  
    outname = fname.GetName()

    if isinstance(dirname,list) or isinstance(histname,list):
      if not isinstance(dirname,list): dirname = [dirname]
      if not isinstance(histname,list): histname = [histname]

      firstHist=True
      for d in dirname: 
        for h in histname: 
          htemp = getHistogram(fname, h, d, postfitmode, allowEmpty, logx)[0]
          if firstHist: 
            histo=htemp.Clone()
            firstHist=False
          else: histo.Add(htemp)
      return [histo,outname]
    
    for key in fname.GetListOfKeys():
        histo = fname.Get(key.GetName())
        dircheck = False
        if dirname == '' : dircheck=True
        elif dirname in key.GetName(): dircheck=True
        if isinstance(histo,ROOT.TH1F) and key.GetName()==histname:
            if logx:
                bin_width = histo.GetBinWidth(1)
                xbins = []
                xbins.append(bin_width - 1)
                axis = histo.GetXaxis()
                for i in range(1,histo.GetNbinsX()+1):
                    xbins.append(axis.GetBinUpEdge(i))
                rethist = ROOT.TH1F(histname,histname,histo.GetNbinsX(),array('d',xbins))
                rethist.SetBinContent(1,histo.GetBinContent(1)*(histo.GetBinWidth(1)-(bin_width - 1))/(histo.GetBinWidth(1)))
                rethist.SetBinError(1,histo.GetBinError(1)*(histo.GetBinWidth(1)-(bin_width - 1))/(histo.GetBinWidth(1)))
                for i in range(2,histo.GetNbinsX()+1):
                    rethist.SetBinContent(i,histo.GetBinContent(i))
                    rethist.SetBinError(i,histo.GetBinError(i))
                histo = rethist
            return [histo,outname]
        elif isinstance(histo,ROOT.TDirectory) and postfitmode in key.GetName() and dircheck:
            return getHistogram(histo,histname, allowEmpty=allowEmpty, logx=logx)
    print 'Failed to find %(postfitmode)s histogram with name %(histname)s in file %(fname)s '%vars()
    if allowEmpty:
        return [ROOT.TH1F('empty', '', 1, 0, 1), outname]
    else:
        return None


def signalComp(leg,plots,colour,stacked):
    return dict([('leg_text',leg),('plot_list',plots),('colour',colour),('in_stack',stacked)])

def backgroundComp(leg,plots,colour):
    return dict([('leg_text',leg),('plot_list',plots),('colour',colour)])

def createAxisHists(n,src,xmin=0,xmax=499):
    result = []
    for i in range(0,n):
        res = src.Clone()
        res.Reset()
        res.SetTitle("")
        res.SetName("axis%(i)d"%vars())
        res.SetAxisRange(xmin,xmax)
        res.SetStats(0)
        result.append(res)
    return result

def PositionedLegendUnrolled(width, height, pos, offset):
    o = offset
    w = width
    h = height
    l = ROOT.gPad.GetLeftMargin()
    t = ROOT.gPad.GetTopMargin()
    b = ROOT.gPad.GetBottomMargin()
    r = ROOT.gPad.GetRightMargin()
    if pos == 1:
        return ROOT.TLegend(l + o, 1 - t - o - h, l + o + w, 1 - t - o, '', 'NBNDC')
    if pos == 2:
        c = l + 0.5 * (1 - l - r)
        return ROOT.TLegend(c - 0.5 * w, 1 - t - o - h, c + 0.5 * w, 1 - t - o, '', 'NBNDC')
    if pos == 3:
        return ROOT.TLegend(1 - r - o - w, 1 - t - o - h, 1 - r - o, 1 - t - o, '', 'NBNDC')
    if pos == 4:
        return ROOT.TLegend(l + o, b + o, l + o + w, b + o + h, '', 'NBNDC')
    if pos == 5:
        c = l + 0.5 * (1 - l - r)
        return ROOT.TLegend(c - 0.5 * w, b + o, c + 0.5 * w, b + o + h, '', 'NBNDC')
    if pos == 6:
        return ROOT.TLegend(1 - r - o - w, b + o, 1 - r - o, b + o + h, '', 'NBNDC')
    if pos == 7:
        return ROOT.TLegend(1 - o - w, 1 - t - o - h, 1 - o, 1 - t - o, '', 'NBNDC')

def DrawTitleUnrolled(pad, text, align, scale=1):
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

    y_off = 1 - t + textOffset * t + 0.01
    if align == 1:
        latex.SetTextAlign(11)
        latex.DrawLatex(l, y_off, text)
    if align == 2:
        latex.SetTextAlign(21)
        latex.DrawLatex(l + (1 - l - r) * 0.5, y_off, text)
    if align == 3:
        latex.SetTextAlign(31)
        latex.DrawLatex(1 - r, y_off, text)
    pad_backup.cd()

    
def parse_arguments():
    parser = argparse.ArgumentParser()
    #Ingredients when output of PostFitShapes is already provided
    parser.add_argument('--file', '-f',
                    help='Input file if shape file has already been created')
    parser.add_argument('--channel',default='',
                    help='Option to specify channel in case it is not obtainable from the shape file name')
    parser.add_argument('--file_dir',default='',
                    help='Name of TDirectory inside shape file')
    parser.add_argument('--mode',default='postfit',
                    help='Prefit or postfit')
    #Blinding options
    parser.add_argument('--manual_blind', action='store_true',
                    default=False,help='Blind data with hand chosen range')
    parser.add_argument('--x_blind_min',default=70,
                    help='Minimum x for manual blinding')
    parser.add_argument('--x_blind_max',default=110,
                    help='Maximum x for manual blinding')
    parser.add_argument('--empty_bin_error',action='store_true',
                    default=False, help='Draw error bars for empty bins')
    #General plotting options
    parser.add_argument('--ratio', default=False,action='store_true',
                    help='Draw ratio plot')
    parser.add_argument('--custom_x_range', action='store_true', 
                    default=False, help='Fix x axis range')
    parser.add_argument('--x_axis_min', default=0.0,
                    help='Fix x axis minimum')
    parser.add_argument('--x_axis_max', default=1000.0, 
                    help='Fix x axis maximum')
    parser.add_argument('--custom_y_range', action='store_true', 
                    default=False, help='Fix y axis range')
    parser.add_argument('--y_axis_min', action='store', default=1., 
                    help='Fix y axis minimum')
    parser.add_argument('--y_axis_max', action='store', default=100.0,
                    help='Fix y axis maximum')
    parser.add_argument('--log_y', action='store_true',
                    help='Use log for y axis')
    parser.add_argument('--log_x', action='store_true',
                    help='Use log for x axis')
    parser.add_argument('--extra_pad', default=0.0, 
                    help='Fraction of extra whitespace at top of plot')
    parser.add_argument('--outname',default='',
                    help='Optional string for start of output filename')
    parser.add_argument('--ratio_range',  default="", 
                    help='y-axis range for ratio plot in format MIN,MAX')
    parser.add_argument('--no_signal', action='store_true',
                    help='Do not draw signal')
    parser.add_argument('--use_asimov', default=False, 
                    action='store_true', help='')
    parser.add_argument('--proper_errors_uniform', default=False, 
                    action='store_true', help='')

    return parser.parse_args()

def main(args):

    fitvars='m_sv_vs_pt_tt'
    fitvars='m_sv'

    era = args.file_dir.split("_")[3]
    if era == "2016":
        lumi = "36.3 fb^{-1} (13 TeV)"
    elif era == "2017":
        lumi = "41.5 fb^{-1} (13 TeV)"
    elif era == "2018":
        lumi = "59.7 fb^{-1} (13 TeV)"
    elif era == "all":
        lumi = "138 fb^{-1} (13 TeV)"

    plot.ModTDRStyle(width=1800, height=700, r=0.4, l=0.16, t=0.12,b=0.15)
    ROOT.TGaxis.SetExponentOffset(-0.06, 0.01, "y")
    # Channel & Category label

    if args.channel == '':
        args.channel = args.file_dir.split("_")[1]
    if args.channel == "tt":
        channel_label = "#tau_{h}#tau_{h}"
    if args.channel == "mt":
        channel_label = "#mu_{}#tau_{h}"
    if args.channel == "et":
        channel_label = "e_{}#tau_{h}"
    if args.channel == "em":
        channel_label = "e_{}#mu_{}"

    bin_number = args.file_dir.split("_")[2]
    #if args.ratio_range=="":
    #  args.ratio_range = "0.7,1.3"
    if bin_number in ["2","35","36","37"]:
        if args.channel=='tt': bin_label = "B-tag"
        if args.channel in ['mt','et']: 
          bin_label = "B-tag"
          if bin_number=="35": bin_label = "B-tag Tight-m_{T}"
          if bin_number=="36": bin_label = "B-tag Loose-m_{T}"
        if args.channel in ['em']:
          bin_label = "B-tag"
          if bin_number=="2": bin_label = "t#bar{t} CR"
          if bin_number=="35": bin_label = "B-tag High d_{#zeta}"
          if bin_number=="36": bin_label = "B-tag Medium d_{#zeta}" 
          if bin_number=="37": bin_label = "B-tag Low d_{#zeta}" 
        plot.ModTDRStyle(r=0.04, l=0.18)
        if args.ratio_range=="":
          args.ratio_range = "0.85,1.15"
    if bin_number in ["32","33","34"]:
        if args.ratio_range=="":
          args.ratio_range = "0.4,1.6"
        if args.channel=='tt': bin_label = "No B-tag"
        if args.channel in ['mt','et']: 
          bin_label = "No B-tag"
          if bin_number=="32": bin_label = "No B-tag Tight-m_{T}"
          if bin_number=="33": bin_label = "No B-tag Loose-m_{T}"
        if args.channel in ['em']:
          bin_label = "No B-tag"
          if bin_number=="32": bin_label = "No B-tag High d_{#zeta}"
          if bin_number=="33": bin_label = "No B-tag Medium d_{#zeta}"
          if bin_number=="34": bin_label = "No B-tag Low d_{#zeta}"

    if bin_number in ["132","232","332","432","133","233","333","433"]:
        if args.ratio_range=="":
          args.ratio_range = "0.7,1.3"
        if args.channel=='tt': bin_label = "No B-tag"
        if args.channel in ['mt','et']:
          bin_label = "No B-tag"
          if bin_number[1:]=="32": bin_label = "No B-tag Tight-m_{T}"
          if bin_number[1:]=="33": bin_label = "No B-tag Loose-m_{T}"
        if args.channel in ['em']:
          bin_label = "No B-tag"
          if bin_number[1:]=="32": bin_label = "No B-tag High d_{#zeta}"
          if bin_number[1:]=="33": bin_label = "No B-tag Medium d_{#zeta}"
          if bin_number[1:]=="34": bin_label = "No B-tag Low d_{#zeta}"

        if bin_number[0] == '1': bin_label+=', p_{T}<50 GeV'
        if bin_number[0] == '2': bin_label+=', 50#leq p_{T}<100 GeV'
        if bin_number[0] == '3': bin_label+=', 100#leq p_{T}<200 GeV'
        if bin_number[0] == '4': bin_label+=', p_{T}#geq 200 GeV'
        plot.ModTDRStyle(r=0.04, l=0.18)

    ## Add bin labels
    bin_labels = {}
    with open("scripts/bin_labels.json") as jsonfile:
          full_bin_labels = json.load(jsonfile)
          if bin_number in ["32","33","34"]:
            bin_labels = full_bin_labels[fitvars]['nobtag']
          else:
            bin_labels= full_bin_labels[fitvars]['btag']

    is2D=False

    print fitvars, bin_number
    if fitvars=='m_sv' or bin_number in ["2","35","36","37","132","232","332","432","133","233","333","433"]:
        x_title = "m_{#tau#tau} (GeV)"
        #x_bins = re.split("\[|\]",bin_labels)[1].split(",")
        #Nxbins = len(x_bins) - 1
    else:
        is2D=True 
        x_title = "Bin number"
        x_bins = re.split("\[|\]",bin_labels)[3].split(",")
        Nxbins = len(x_bins) - 1

    if is2D:
        y_bin_var = re.split(",|\[|\]", bin_labels)[0]
        y_bin_labels = re.split("\[|\]",bin_labels)[1].split(",")
    else:
        y_bin_var = ""
        y_bin_labels = ""

    is2D = False # we now split 2D histograms into 1D

    file_dir = args.file_dir
    mode = args.mode
    manual_blind = args.manual_blind
    x_blind_min = args.x_blind_min
    x_blind_max = args.x_blind_max
    empty_bin_error = args.empty_bin_error
    extra_pad = float(args.extra_pad)
    custom_x_range = args.custom_x_range
    custom_y_range = args.custom_y_range
    x_axis_min = float(args.x_axis_min)
    x_axis_max = float(args.x_axis_max)
    y_axis_min = float(args.y_axis_min)
    y_axis_max = float(args.y_axis_max)
    log_y=args.log_y
    log_x=args.log_x
    if(args.outname != ''):
        outname=args.outname + '_'
    else:
        outname=''
    
    
    if args.file:
        print "Providing shape file: ", args.file, ", with specified subdir name: ", file_dir
        shape_file=args.file
        shape_file_name=args.file
    
    histo_file = ROOT.TFile(shape_file)
    
    #Store plotting information for different backgrounds 
    background_schemes = {
        'mt':[
                backgroundComp("H(125 GeV)#rightarrow#tau#tau",["qqH125","bbH125","ggH125"],ROOT.TColor.GetColor(51,51,230)),
                backgroundComp("Diboson",["VVL"],ROOT.TColor.GetColor("#6F2D35")),
                backgroundComp("t#bar{t}",["TTL"],ROOT.TColor.GetColor("#9999cc")),
                backgroundComp("Z#rightarrowll",["ZL"],ROOT.TColor.GetColor("#4496c8")),
                backgroundComp("Jet#rightarrow#tau_{h}",["jetFakes"],ROOT.TColor.GetColor(192, 232, 100)),
                backgroundComp("#mu#rightarrow#tau embedding",["EMB"],ROOT.TColor.GetColor("#ffcc66")),
                ],
        'et':[
                backgroundComp("H(125 GeV)#rightarrow#tau#tau",["qqH125","bbH125","ggH125"],ROOT.TColor.GetColor(51,51,230)),
                backgroundComp("Diboson",["VVL"],ROOT.TColor.GetColor("#6F2D35")),
                backgroundComp("t#bar{t}",["TTL"],ROOT.TColor.GetColor("#9999cc")),
                backgroundComp("Z#rightarrowll",["ZL"],ROOT.TColor.GetColor("#4496c8")),
                backgroundComp("Jet#rightarrow#tau_{h}",["jetFakes"],ROOT.TColor.GetColor(192, 232, 100)),
                backgroundComp("#mu#rightarrow#tau embedding",["EMB"],ROOT.TColor.GetColor("#ffcc66")),
                ],
        'tt':[
                backgroundComp("H(125 GeV)#rightarrow#tau#tau",["qqH125","bbH125","ggH125"],ROOT.TColor.GetColor(51,51,230)),
                backgroundComp("Diboson",["VVL"],ROOT.TColor.GetColor("#6F2D35")),
                backgroundComp("t#bar{t}",["TTL"],ROOT.TColor.GetColor("#9999cc")),
                backgroundComp("Z#rightarrowll",["ZL"],ROOT.TColor.GetColor("#4496c8")),
                backgroundComp("Jet#rightarrow#tau_{h}",["jetFakes","wFakes"],ROOT.TColor.GetColor(192, 232, 100)),
                backgroundComp("#mu#rightarrow#tau embedding",["EMB"],ROOT.TColor.GetColor("#ffcc66")),
                ],


        'em':[
                backgroundComp("H(125 GeV)#rightarrow#tau#tau",["qqH125","bbH125","ggH125"],ROOT.TColor.GetColor(51,51,230)),
                backgroundComp("QCD", ["QCD"], ROOT.TColor.GetColor("#ffccff")),
                backgroundComp("Diboson",["VVL"],ROOT.TColor.GetColor("#6F2D35")),
                backgroundComp("W+jets",["W"],ROOT.TColor.GetColor(222, 90, 106)),
                backgroundComp("t#bar{t}",["TTL"],ROOT.TColor.GetColor("#9999cc")),
                backgroundComp("Z#rightarrowll",["ZL"],ROOT.TColor.GetColor("#4496c8")),
                backgroundComp("#mu#rightarrow#tau embedding",["EMB"],ROOT.TColor.GetColor("#ffcc66")),
                ],

        }

    #Extract relevent histograms from shape file
    sighists = []

    #signal_names = 'TotalSig'
    signal_names = ['ggh_t','ggh_i','ggh_b']
    if args.no_signal:  signal_names = []

    file_dir_list = []
    file_dir_list = [file_dir]

    if not args.no_signal: 
      [sighist,binname] = getHistogram(histo_file,signal_names, file_dir_list, mode, args.no_signal, log_x)
      sighists.append(sighist)
    bkghist = getHistogram(histo_file,'TotalBkg',file_dir, mode, logx=log_x)[0]
    sbhist = bkghist.Clone()
    # can use this one for showing ggX as well as qqX
    sbhist_alt = bkghist.Clone()
    
    if not args.use_asimov:
        total_datahist = getHistogram(histo_file,"data_obs",file_dir, mode, logx=log_x)[0]
    else:
        total_datahist = getHistogram(histo_file,"TotalProcs",file_dir, mode, logx=log_x)[0].Clone()
        for bin_ in range(1,total_datahist.GetNbinsX()+1):
            content = total_datahist.GetBinContent(bin_)
            total_datahist.SetBinError(bin_, np.sqrt(content))


    blind_datahist = total_datahist.Clone()
    total_datahist.SetMarkerStyle(20)
    blind_datahist.SetMarkerStyle(20)
    blind_datahist.SetLineColor(1)

    total_bkg = getHistogram(histo_file,"TotalBkg",file_dir, mode, logx=log_x)[0].Clone()
    azimov_datahist = blind_datahist.Clone()
    for i in range(0, azimov_datahist.GetNbinsX()+1):
      azimov_datahist.SetBinContent(i,-0.1)
      azimov_datahist.SetBinError(i,0)

    azimov_datahist.SetLineColor(ROOT.kRed)
    azimov_datahist.SetMarkerColor(ROOT.kRed)

    #Blinding by hand using requested range, set to 70-110 by default
    # for 0jet category
    if not is2D and manual_blind:
        for i in range(0,total_datahist.GetNbinsX()):
            low_edge = total_datahist.GetBinLowEdge(i+1)
            high_edge = low_edge+total_datahist.GetBinWidth(i+1)
            if ((low_edge > float(x_blind_min) and low_edge < float(x_blind_max)) 
                    or (high_edge > float(x_blind_min) and high_edge<float(x_blind_max))):
                blind_datahist.SetBinContent(i+1, -0.1)
                blind_datahist.SetBinError(i+1,0)
                c = total_bkg.GetBinContent(i+1)
                azimov_datahist.SetBinContent(i+1,c)
                azimov_datahist.SetBinError(i+1,c**.5)
    # for boosted category:
    if is2D and manual_blind:
        x_blind_ind = [ind for ind, x in enumerate(x_bins) if 120 >= int(x) >= 70]
        x_blind_ind1 = []
        
        dummy_list = [int(x) for x in np.linspace(Nxbins,Nxbins*Nxbins,Nxbins)]
        for i in range(0,total_datahist.GetNbinsX()):
            if i in dummy_list:
                x_blind_ind1 = [x+i for x in x_blind_ind]
            if i in x_blind_ind or i in x_blind_ind1:
                blind_datahist.SetBinContent(i+1,-0.1)
                blind_datahist.SetBinError(i+1,0)
                c = total_bkg.GetBinContent(i+1)
                azimov_datahist.SetBinContent(i+1,c)
                azimov_datahist.SetBinError(i+1,c**.5)


    #Set bin errors for empty bins if required:
    if empty_bin_error:
        for i in range (1,blind_datahist.GetNbinsX()+1):
            if blind_datahist.GetBinContent(i) == 0:
                blind_datahist.SetBinError(i,1.8)
    #Set uniform bin errors properly for Content < 10 bins
    if args.proper_errors_uniform:
        proper_errs_dict = {
                0: 1.29, 1: 2.38, 2: 3.51, 3: 4.20, 4: 4.44, 5: 5.06,
                6: 5.46, 7: 6.05, 8: 6.02, 9: 6.46 
                }
        for i in range (1,blind_datahist.GetNbinsX()+1):
            if blind_datahist.GetBinContent(i) < 9.5 and blind_datahist.GetBinContent(i) >= 0:
                new_err = proper_errs_dict[round(blind_datahist.GetBinContent(i))]
                blind_datahist.SetBinError(i, new_err)

    #Normalise by bin width 
    scale=1.0
    if is2D: scale=1./10. # if we have an unrolled plot then we need to account for the bins being in 10 GeV units
    bkghist.Scale(scale,"width")
    sbhist.Scale(scale,"width")

    for shist in sighists:
        shist.Scale(5.4*scale,"width") # can scale up signals here if desired
        sbhist.Add(shist)
    #sbhist.Scale(scale,"width")
    #sbhist_alt.Scale(scale,"width")
    for shist in sighists:
        shist.Scale(scale,"width")
    blind_datahist.Scale(scale,"width")
    azimov_datahist.Scale(scale,"width")

    blind_datagraph = ROOT.TGraphAsymmErrors(blind_datahist)
    azimov_datagraph = ROOT.TGraphAsymmErrors(azimov_datahist)

    channel = args.channel
    
    #Create stacked plot for the backgrounds
    bkg_histos = []
    for i,t in enumerate(background_schemes[channel]):
        plots = t['plot_list']
        isHist = False
        h = ROOT.TH1F()
        for j,k in enumerate(plots):
            if h.GetEntries()==0 and getHistogram(histo_file,k, file_dir,mode,False,logx=log_x) is not None:
                isHist = True
                h = getHistogram(histo_file,k, file_dir,mode, logx=log_x)[0]
                h.SetName(k)
            else:
                if getHistogram(histo_file,k, file_dir,mode, False, logx=log_x) is not None:
                    isHist = True
                    h.Add(getHistogram(histo_file,k, file_dir,mode,logx=log_x)[0])
        h.SetFillColor(t['colour'])
        h.SetLineColor(ROOT.kBlack)
        h.SetMarkerSize(0)
        
        h.Scale(scale,"width")
        if isHist:
            bkg_histos.append(h)
    
    stack = ROOT.THStack("hs","")
    for hists in bkg_histos:
        stack.Add(hists)
    
    #Setup style related things
    c2 = ROOT.TCanvas()
    c2.cd()
    
    if args.ratio:
        pads=plot.TwoPadSplit(0.35,0.01,0.01)
    else:
        pads=plot.OnePad()

    for p in pads: p.SetFrameLineWidth(1)

    pads[0].cd()
    if(log_y):
        pads[0].SetLogy(1)
    if(log_x):
        pads[0].SetLogx(1)
    
    if custom_x_range:
            if x_axis_max > bkghist.GetXaxis().GetXmax(): x_axis_max = bkghist.GetXaxis().GetXmax()
    if args.ratio:
        if(log_x): 
            pads[1].SetLogx(1)
        axish = createAxisHists(2,bkghist,bkghist.GetXaxis().GetXmin(),bkghist.GetXaxis().GetXmax()-0.01)
        axish[1].GetXaxis().SetTitle(x_title)
        if is2D:
            axish[1].GetXaxis().SetLabelSize(0.03)
        axish[1].GetYaxis().SetNdivisions(4)
        axish[1].GetYaxis().SetLabelSize(0.033)
        axish[1].GetXaxis().SetLabelSize(0.033)
        if is2D:
            axish[1].GetXaxis().SetNdivisions(bkghist.GetNbinsX()/Nxbins,Nxbins,0,False)
        axish[0].GetYaxis().SetTitleSize(0.048)
        axish[0].GetYaxis().SetLabelSize(0.033)
        axish[0].GetYaxis().SetTitleOffset(0.5)
        axish[0].GetXaxis().SetTitleSize(0)
        axish[0].GetXaxis().SetLabelSize(0)
        if is2D:
            axish[0].GetXaxis().SetNdivisions(bkghist.GetNbinsX()/Nxbins,Nxbins,0,False)
        axish[0].GetYaxis().SetTitleSize(axish[1].GetXaxis().GetTitleSize())
        axish[0].GetXaxis().SetRangeUser(x_axis_min,bkghist.GetXaxis().GetXmax()-0.01)
        axish[1].GetXaxis().SetRangeUser(x_axis_min,bkghist.GetXaxis().GetXmax()-0.01)
        axish[0].GetXaxis().SetMoreLogLabels()
        axish[0].GetXaxis().SetNoExponent()
        axish[1].GetXaxis().SetMoreLogLabels()
        axish[1].GetXaxis().SetNoExponent()
        axish[1].GetXaxis().SetTitleOffset(0.85)

        if not is2D:
            axish[0].GetXaxis().SetRangeUser(0.,bkghist.GetXaxis().GetXmax()-0.01)
            axish[1].GetXaxis().SetRangeUser(0.,bkghist.GetXaxis().GetXmax()-0.01)
 
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min,x_axis_max-0.01)
            axish[1].GetXaxis().SetRangeUser(x_axis_min,x_axis_max-0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min,y_axis_max)
            axish[1].GetYaxis().SetRangeUser(y_axis_min,y_axis_max)
    else:
        axish = createAxisHists(1,bkghist,bkghist.GetXaxis().GetXmin(),bkghist.GetXaxis().GetXmax()-0.01)
        if custom_x_range:
            axish[0].GetXaxis().SetRangeUser(x_axis_min,x_axis_max-0.01)
        if custom_y_range:
            axish[0].GetYaxis().SetRangeUser(y_axis_min,y_axis_max)
    axish[0].GetYaxis().SetTitleOffset(0.5)
    axish[1].GetYaxis().SetTitleOffset(0.5)
    axish[0].GetYaxis().SetTitle("dN/dm_{#tau#tau} (1/GeV)")
    if not is2D:
        axish[0].GetYaxis().SetTitle("dN/dm_{#tau#tau} (1/GeV)")
        axish[0].GetYaxis().SetTitleOffset(1.5)
        axish[1].GetYaxis().SetTitleOffset(1.5)

    axish[0].GetXaxis().SetTitle(x_title)
        
    if not custom_y_range: axish[0].SetMaximum(extra_pad*bkghist.GetMaximum())
    if not custom_y_range: 
        if(log_y): 
            ymin = 0.1
            axish[0].SetMinimum(ymin)
        else: 
            axish[0].SetMinimum(0)
    
    hist_indices = [0] 
    for i in hist_indices:
        pads[i].cd()
        axish[i].Draw("AXIS")

        bkghist.SetMarkerSize(0)
        bkghist.SetFillColor(2001)
        bkghist.SetLineColor(1)
        bkghist.SetLineWidth(1)
        bkghist.SetFillColor(plot.CreateTransparentColor(12,0.4))
        bkghist.SetLineColor(0)    

        stack.Draw("histsame")
        bkghist.Draw("e2same")
        #Add signal, either model dependent or independent
        if not args.no_signal:
            sighist.SetLineColor(ROOT.kRed)
            sighist.SetLineWidth(2)
            if not is2D:
                sighist.SetLineWidth(3)
            for shist in sighists:
                for j in range(1,shist.GetNbinsX()+1):
                    entry = shist.GetBinContent(j)
                    if entry < axish[0].GetMinimum():
                        shist.SetBinContent(j,axish[0].GetMinimum()*1.00001)
                shist.Draw("histsame][") # removing vertical lines at the borders of the pad; possible with the trick above
        blind_datagraph_extra = blind_datagraph.Clone()
        blind_datagraph_extra.Draw("P Z 0 same")
        blind_datagraph.SetMarkerSize(0.)
        blind_datagraph.Draw("P Z 0 same")

        azimov_datagraph_extra = azimov_datagraph.Clone()
        azimov_datagraph_extra.Draw("P Z 0 same")
        azimov_datagraph.SetMarkerSize(0.)
        azimov_datagraph.Draw("P Z 0 same")

        axish[i].Draw("axissame")
    
    pads[0].cd()
    pads[0].SetTicks(1)
    pads[1].SetTicks(1)
    #Setup legend
    if not is2D:
        legend = plot.PositionedLegend(0.35,0.30,3,0.03)
        legend.SetTextSize(0.025)
    else:
        legend = PositionedLegendUnrolled(0.13,0.5,7,0.02)
        legend.SetTextSize(0.035) 
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    
    legend.AddEntry(total_datahist,"Observed","PEl")
    #Drawn on legend in reverse order looks better
    bkg_histos.reverse()

    background_schemes[channel].reverse()
    leg_hists = [None]*len(bkg_histos)
    for legi,hists in enumerate(bkg_histos):
        legend.AddEntry(hists,background_schemes[channel][legi]['leg_text'],"f")
    #legend.AddEntry(bkghist,"Background uncertainty","f")
    bkghist.SetLineWidth(0)
    legend.AddEntry(bkghist,"Bkg. Uncertainty","f")
    if is2D:
        if not mode == 'prefit':
          if not args.no_signal: legend.AddEntry(sighist,"qqH(95 GeV) @ 1 pb"%vars(),"l")
        else:
          if not args.no_signal: legend.AddEntry(sighist,"qqH(95 GeV) @ 1 pb"%vars(),"l")

    else:
        if not args.no_signal: legend.AddEntry(sighist,"qqH(95 GeV) @ 1 pb"%vars(),"l")
    legend.Draw("same")

    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(ROOT.kBlack)
    latex2.SetTextFont(42)
    if not is2D:
        latex2.SetTextSize(0.032)
        latex2.DrawLatex(0.21,0.89,"{}, {}".format(bin_label, channel_label))
    else:
        latex2.SetTextAlign(23)
        latex2.SetTextSize(0.05)
        latex2.DrawLatex(0.46,0.955,"{}, {}".format(bin_label, channel_label))

    #CMS and lumi labels
    plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), extra_pad if extra_pad>0 else 0.15)
    extra=''
    extra='Preliminary'
    if 'htt_tt_2018_6_' in file_dir: extra=''
    if not is2D:
        DrawCMSLogo(pads[0], 'CMS', extra, 0, 0.07, -0.0, 2.0, '', 0.85, relExtraDX=0.05)
        plot.DrawTitle(pads[0], lumi, 3, textSize=0.6)
    else:
        DrawCMSLogo(pads[0], 'CMS', extra, 0, 0.07, -0.0, 2.0, '', 0.6, relExtraDX=0.005)
        DrawTitleUnrolled(pads[0], lumi, 3, scale=0.6)
    
    #Add ratio plot if required
    axish[1].GetYaxis().SetTitle("Obs./Exp.")
    if args.ratio:
        ratio_bkghist = plot.MakeRatioHist(bkghist,bkghist,True,False)
        sbhist.SetLineColor(ROOT.kRed)
        sbhist_alt.SetLineColor(ROOT.kBlue)
        sbhist.SetLineWidth(2)
        sbhist_alt.SetLineWidth(2)
        if not is2D:
            sbhist.SetLineWidth(3)
            sbhist_alt.SetLineWidth(3)

        bkghist_errors = bkghist.Clone()
        #for i in range(1,bkghist_errors.GetNbinsX()+1): bkghist_errors.SetBinContent(i,bkghist.GetBinError(i))
        for i in range(1,bkghist_errors.GetNbinsX()+1): bkghist_errors.SetBinContent(i,1.)
        ratio_sighist = plot.MakeRatioHist(sbhist,bkghist,True,False)
        ratio_datahist_ = plot.MakeRatioHist(blind_datahist,bkghist,True,False)
        azimov_ratio_datahist_ = plot.MakeRatioHist(azimov_datahist,bkghist,True,False)

        ratio_datahist = ROOT.TGraphAsymmErrors(ratio_datahist_)
        azimov_ratio_datahist = ROOT.TGraphAsymmErrors(azimov_ratio_datahist_)

        pads[1].cd()
        pads[1].SetGrid(0,1)
        axish[1].Draw("axis")
        axish[1].SetMinimum(float(args.ratio_range.split(',')[0]))
        axish[1].SetMaximum(float(args.ratio_range.split(',')[1]))
        ratio_bkghist.SetMarkerSize(0)

        ratio_bkghist.Draw("e2same")
        ratio_datahist.Draw("P Z 0 same")
        ratio_sighist.Draw('histsame')
        if args.manual_blind: azimov_ratio_datahist.Draw("P Z 0 same")
        pads[1].RedrawAxis("G")
        if is2D:
            rlegend = ROOT.TLegend(0.85, 0.27, 0.98, 0.16, '', 'NBNDC')
            rlegend.SetTextFont(42)
            rlegend.SetTextSize(0.035)
            rlegend.SetFillStyle(0)
            rlegend.AddEntry(ratio_datahist,"Obs./Exp.","PE")
            rlegend.AddEntry(""," ","")
            #rlegend.AddEntry(ratio_sighist,"(Sig.+Bkg.)/Bkg.","L")
        #rlegend.Draw("same")
        # Draw extra axis for explanation (use "N" for no optimisation)
        if is2D:
            extra_axis = ROOT.TGaxis(0,-0.03,30.,-0.03,0.,300.,403,"NS")
            extra_axis.SetLabelSize(0.03)
            extra_axis.SetLabelFont(42)
            extra_axis.SetMaxDigits(3)
            extra_axis.SetTitle("m_{#tau#tau} (GeV)")
            extra_axis.SetTitleFont(42)
            extra_axis.SetTitleSize(0.035)
            extra_axis.SetTitleOffset(1.1)
            extra_axis.SetTickSize(0.08)
            extra_axis.Draw()

    pads[0].cd()
    pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    ## Add lines after every Nxbins 
    #line = ROOT.TLine()
    #line.SetLineWidth(2)
    #line.SetLineStyle(3)
    #line.SetLineColor(1)
    #x = bkghist.GetNbinsX()/Nxbins
    #x_min=bkghist.GetXaxis().GetBinLowEdge(1)
    #x_max=bkghist.GetXaxis().GetBinLowEdge(bkghist.GetNbinsX()+1)
    #x_range=x_max-x_min

    ## for now we have to hard code x=4 otherwise we have issues when the auto rebinning has been used
    #x=4
    #if is2D:
    #    for l in range(1,x):
    #        pads[0].cd()
    #        ymax = axish[0].GetMaximum()
    #        ymin = axish[0].GetMinimum()
    #        line.DrawLine(l*x_range/x,ymin,l*x_range/x,ymax)
    #        if args.ratio:
    #            pads[1].cd()
    #            ymax = axish[1].GetMaximum()
    #            ymin = axish[1].GetMinimum()
    #            line.DrawLine(l*x_range/x,ymin,l*x_range/x,ymax)
    
    ## Add bin labels between lines
    pads[0].cd()
    latex_bin = ROOT.TLatex()
    latex_bin.SetNDC()
    latex_bin.SetTextAngle(0)
    latex_bin.SetTextColor(ROOT.kBlack)
    latex_bin.SetTextFont(42)
    latex_bin.SetTextSize(0.035)
    if len(y_bin_labels) > 6 or ('_2_' in file_dir and len(y_bin_labels)>5):
        latex_bin.SetTextSize(0.027)

    for i in range(0, len(y_bin_labels)):
        if i < len(y_bin_labels)-1:
            y_bin_label = "{} #leq {} < {} {}".format(y_bin_labels[i],y_bin_var,y_bin_labels[i+1],"GeV")
            l = ROOT.gPad.GetLeftMargin()
            r = ROOT.gPad.GetRightMargin()
            #xshift = ((1-r-l)/(axish[0].GetNbinsX()/Nxbins))*i + l
            # have to hard code 4 here other wise we have problems when the autorebinning has been used
            xshift = ((1-r-l)/(4))*i + l
            latex_bin.DrawLatex(xshift+0.02,0.82,y_bin_label)
        else:
            #xshift = ((1-r-l)/(axish[0].GetNbinsX()/Nxbins))*i + l
            # have to hard code 4 here other wise we have problems when the autorebinning has been used
            xshift = ((1-r-l)/(4))*i + l
            y_bin_label = "{} > {} {}".format(y_bin_var,y_bin_labels[i],"GeV")
            latex_bin.DrawLatex(xshift+0.02,0.82,y_bin_label)
    
    #Save as png and pdf with some semi sensible filename
    shape_file_name = shape_file_name.replace(".root","_%(mode)s"%vars())
    shape_file_name = shape_file_name.replace("_shapes","")
    outname += shape_file_name+"_"+file_dir.strip("htt").strip("_")
    if(log_x): 
        outname+="_logx"
    c2.SaveAs("plots/%(outname)s.pdf"%vars())

    del c2
    histo_file.Close()

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
    

