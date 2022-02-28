import ROOT
ROOT.gROOT.SetBatch(True)
import CombineHarvester.CombineTools.plotting as plot

# before you can use this plotting script you need to run the postfit shapes like:
# for x in plot_50to100 plot_100to200 plot_0to50 plot_GT200; do PostFitShapesFromWorkspace -w model_independent_limits/Dec01_2D_all_all_bsm-model-indep/combined/plot_${x}/ws.root -d model_independent_limits/Dec01_2D_all_all_bsm-model-indep/combined/plot_${x}/combined.txt.cmb --fitresult model_independent_limits/Nov19_2D_all_all_bsm-model-indep/combined/cmb/multidimfit.ggH.m100.bestfit.robustHesse.root:fit_mdf -o shapes_prop_plot_postfit_plot_${x}_v5.root --skip-prefit=true   --mass 100 --total-shapes=true --postfit; done 

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

COL_STORE = []

c1 = ROOT.TCanvas()

cats = ['_0to50','_50to100','_100to200','_GT200','']

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



def SubtractData(h_b, h_d):
  out = h_d.Clone()
  out.SetName(out.GetName()+'_sub')
  for i in range(1,h_b.GetNbinsX()+1):
    c = h_d.GetBinContent(i) - h_b.GetBinContent(i)
    out.SetBinContent(i,c)
  return out

def DivideData(h_b, h_d):
  out = h_d.Clone()
  out.SetName(out.GetName()+'_ratio')
  for i in range(1,h_b.GetNbinsX()+1):
    c = h_d.GetBinContent(i)/h_b.GetBinContent(i)
    e = h_d.GetBinError(i)/h_b.GetBinContent(i)
    out.SetBinContent(i,c)
    out.SetBinError(i,e)
  return out

def backgroundComp(leg,plots,colour):
  return dict([('leg_text',leg),('plot_list',plots),('colour',colour)])

def TwoPadSplit(split_point, gap_low, gap_high):
    upper = ROOT.TPad('upper', 'upper', 0., 0., 1., 1.)
    upper.SetBottomMargin(split_point + gap_high)
    upper.SetFillStyle(4000)
    upper.SetTicks(1)
    upper.Draw()
    lower = ROOT.TPad('lower', 'lower', 0., 0., 1., 1.)
    lower.SetTopMargin(1 - split_point + gap_low)
    lower.SetFillStyle(4000)
    lower.Draw()
    upper.cd()
    result = [upper, lower]
    return result

def CreateTransparentColor(color, alpha):
    adapt = ROOT.gROOT.GetColor(color)
    new_idx = ROOT.gROOT.GetListOfColors().GetLast() + 1
    trans = ROOT.TColor(
        new_idx, adapt.GetRed(), adapt.GetGreen(), adapt.GetBlue(), '', alpha)
    COL_STORE.append(trans)
    trans.SetName('userColor%i' % new_idx)
    return new_idx

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

def SetTDRStyle():
    """Sets the PubComm recommended style

    Just a copy of <http://ghm.web.cern.ch/ghm/plots/MacroExample/tdrstyle.C>
    @sa ModTDRStyle() to use this style with some additional customisation.
    """
    # For the canvas:
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)  # Height of canvas
    ROOT.gStyle.SetCanvasDefW(600)  # Width of canvas
    ROOT.gStyle.SetCanvasDefX(0)    # POsition on screen
    ROOT.gStyle.SetCanvasDefY(0)

    # For the Pad:
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(ROOT.kWhite)
    ROOT.gStyle.SetPadGridX(False)
    ROOT.gStyle.SetPadGridY(False)
    ROOT.gStyle.SetGridColor(0)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)

    # For the frame:
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetFrameBorderSize(1)
    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetFrameFillStyle(0)
    ROOT.gStyle.SetFrameLineColor(1)
    ROOT.gStyle.SetFrameLineStyle(1)
    ROOT.gStyle.SetFrameLineWidth(1)

    # For the histo:
    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)

    ROOT.gStyle.SetEndErrorSize(2)
    # ROOT.gStyle.SetErrorMarker(20)
    # ROOT.gStyle.SetErrorX(0.)

    ROOT.gStyle.SetMarkerStyle(20)

    # For the fit/function:
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetFitFormat('5.4g')
    ROOT.gStyle.SetFuncColor(2)
    ROOT.gStyle.SetFuncStyle(1)
    ROOT.gStyle.SetFuncWidth(1)

    # For the date:
    ROOT.gStyle.SetOptDate(0)

    # For the statistics box:
    ROOT.gStyle.SetOptFile(0)
    ROOT.gStyle.SetOptStat(0)
    # To display the mean and RMS:   SetOptStat('mr')
    ROOT.gStyle.SetStatColor(ROOT.kWhite)
    ROOT.gStyle.SetStatFont(42)
    ROOT.gStyle.SetStatFontSize(0.025)
    ROOT.gStyle.SetStatTextColor(1)
    ROOT.gStyle.SetStatFormat('6.4g')
    ROOT.gStyle.SetStatBorderSize(1)
    ROOT.gStyle.SetStatH(0.1)
    ROOT.gStyle.SetStatW(0.15)

    # Margins:
    ROOT.gStyle.SetPadTopMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.02)

    # For the Global title:
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    # For the axis titles:
    ROOT.gStyle.SetTitleColor(1, 'XYZ')
    ROOT.gStyle.SetTitleFont(42, 'XYZ')
    ROOT.gStyle.SetTitleSize(0.06, 'XYZ')
    ROOT.gStyle.SetTitleXOffset(0.9)
    ROOT.gStyle.SetTitleYOffset(1.25)

    # For the axis labels:

    ROOT.gStyle.SetLabelColor(1, 'XYZ')
    ROOT.gStyle.SetLabelFont(42, 'XYZ')
    ROOT.gStyle.SetLabelOffset(0.007, 'XYZ')
    ROOT.gStyle.SetLabelSize(0.05, 'XYZ')

    # For the axis:

    ROOT.gStyle.SetAxisColor(1, 'XYZ')
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, 'XYZ')
    ROOT.gStyle.SetNdivisions(510, 'XYZ')
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    # Change for log plots:
    ROOT.gStyle.SetOptLogx(0)
    ROOT.gStyle.SetOptLogy(0)
    ROOT.gStyle.SetOptLogz(0)

    # Postscript options:
    ROOT.gStyle.SetPaperSize(20., 20.)

    ROOT.gStyle.SetHatchesLineWidth(5)
    ROOT.gStyle.SetHatchesSpacing(0.05)


def ModTDRStyle(width=600, height=600, t=0.06, b=0.12, l=0.16, r=0.04):
    """Modified version of the tdrStyle

    Args:
        width (int): Canvas width in pixels
        height (int): Canvas height in pixels
        t (float): Pad top margin [0-1]
        b (float): Pad bottom margin [0-1]
        l (float): Pad left margin [0-1]
        r (float): Pad right margin [0-1]
    """
    SetTDRStyle()

    # Set the default canvas width and height in pixels
    ROOT.gStyle.SetCanvasDefW(width)
    ROOT.gStyle.SetCanvasDefH(height)

    # Set the default margins. These are given as fractions of the pad height
    # for `Top` and `Bottom` and the pad width for `Left` and `Right`. But we
    # want to specify all of these as fractions of the shortest length.
    def_w = float(ROOT.gStyle.GetCanvasDefW())
    def_h = float(ROOT.gStyle.GetCanvasDefH())

    scale_h = (def_w / def_h) if (def_h > def_w) else 1.
    scale_w = (def_h / def_w) if (def_w > def_h) else 1.

    def_min = def_h if (def_h < def_w) else def_w

    ROOT.gStyle.SetPadTopMargin(t * scale_h)
    # default 0.05
    ROOT.gStyle.SetPadBottomMargin(b * scale_h)
    # default 0.13
    ROOT.gStyle.SetPadLeftMargin(l * scale_w)
    # default 0.16
    ROOT.gStyle.SetPadRightMargin(r * scale_w)
    # default 0.02
    # But note the new CMS style sets these:
    # 0.08, 0.12, 0.12, 0.04

    # Set number of axis tick divisions
    ROOT.gStyle.SetNdivisions(510, 'XYZ')  # default 510

    # Some marker properties not set in the default tdr style
    ROOT.gStyle.SetMarkerColor(ROOT.kBlack)
    ROOT.gStyle.SetMarkerSize(1.0)

    ROOT.gStyle.SetLabelOffset(0.007, 'YZ')
    # This is an adhoc adjustment to scale the x-axis label
    # offset when we stretch plot vertically
    # Will also need to increase if first x-axis label has more than one digit
    ROOT.gStyle.SetLabelOffset(0.005 * (3. - 2. / scale_h), 'X')

    # In this next part we do a slightly involved calculation to set the axis
    # title offsets, depending on the values of the TPad dimensions and
    # margins. This is to try and ensure that regardless of how these pad
    # values are set, the axis titles will be located towards the edges of the
    # canvas and not get pushed off the edge - which can often happen if a
    # fixed value is used.
    title_size = 0.05
    title_px = title_size * def_min
    label_size = 0.04
    ROOT.gStyle.SetTitleSize(title_size, 'XYZ')
    ROOT.gStyle.SetLabelSize(label_size, 'XYZ')

    ROOT.gStyle.SetTitleXOffset(0.5 * scale_h *
                             (1.2 * (def_h * b * scale_h - 0.6 * title_px)) /
                             title_px)
    ROOT.gStyle.SetTitleYOffset(0.5 * scale_w *
                             (1.2 * (def_w * l * scale_w - 0.6 * title_px)) /
                             title_px)

    # Only draw ticks where we have an axis
    ROOT.gStyle.SetPadTickX(0)
    ROOT.gStyle.SetPadTickY(0)
    ROOT.gStyle.SetTickLength(0.02, 'XYZ')

    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetFillColor(0)

    ROOT.gROOT.ForceStyle()

def HTTMassPlot(cats=[432],
            infile=None, infile2=None, prefit_plot=None
            ):

    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.TH1.AddDirectory(False)

    ModTDRStyle(r=0.04, l=0.14,height=700)

    sig_schemes = ("ggH(100 GeV)#rightarrow#tau#tau @ 5.7 pb", ["TotalSig"])

    background_schemes = [
                backgroundComp("H(125 GeV)",["qqH125","bbH125","ggH125","ggHWW125", "qqHWW125", "WHWW125", "ZHWW125"],ROOT.TColor.GetColor(51,51,230)),
                backgroundComp("Diboson",["VVL"],ROOT.TColor.GetColor("#6F2D35")),
                backgroundComp("t#bar{t}",["TTL"],ROOT.TColor.GetColor("#9999cc")),
                backgroundComp("Z#rightarrowll",["ZL"],ROOT.TColor.GetColor("#4496c8")),
                backgroundComp("Jet#rightarrow#tau",["jetFakes","wFakes","QCD","W"],ROOT.TColor.GetColor(192, 232, 100)),
                backgroundComp("#mu#rightarrow#tau embedding",["EMB"],ROOT.TColor.GetColor("#ffcc66")),]

    total_datahist = infile.Get('postfit/data_obs').Clone()
    total_bkg = infile.Get('postfit/TotalBkg').Clone()
    total_sig = infile.Get('postfit/TotalSig').Clone()
    if infile2: bkgonly = infile2.Get('postfit/TotalBkg').Clone()
    if prefit_plot: prefit = prefit_plot.Get('prefit/TotalBkg').Clone()

    total_datahist.Scale(1.0,"width")
    total_bkg.Scale(1.0,"width")
    total_sig.Scale(1.0,"width")

    #Create stacked plot for the backgrounds
    bkg_histos = []
    

    for i,t in enumerate(background_schemes):
        plots = t['plot_list']
        h = ROOT.TH1F()
          
        for j,k in enumerate(plots):
          for chan in ['et','mt','tt','em']:
            if chan == 'em':
              cats_ = []
              for c in cats:
                cats_.append(c)
                cats_.append(c+1)
            else: cats_ = cats
            for c in cats_:
              for y in [2016,2017,2018]:
                dirname = 'htt_%(chan)s_%(c)s_%(y)s_postfit' % vars() 
                if not (isinstance(infile.Get(dirname+'/'+k),ROOT.TH1D) or isinstance(infile.Get(dirname+'/'+k),ROOT.TH1F)): continue
 
                if h.GetEntries()==0:
                    h = infile.Get(dirname+'/'+k).Clone()
                    h.SetName(k)
                else:
                    h.Add(infile.Get(dirname+'/'+k).Clone())


        h.SetFillColor(t['colour'])
        h.SetLineColor(ROOT.kBlack)
        h.SetMarkerSize(0)
        h.Scale(1.0,"width")
        if h.GetName() == '': continue
        bkg_histos.append(h)

    stack = ROOT.THStack("hs","")
    for hists in bkg_histos:
      stack.Add(hists.Clone())

    total_sig.SetFillColor(0)
    total_sig.SetLineColor(ROOT.kRed)
    total_sig.SetLineWidth(2)
    stack.Add(total_sig)

    c1 = ROOT.TCanvas()
    c1.cd()

    pads=TwoPadSplit(0.39,0.01,0.01)
    pads[0].cd()

    x_title = 'm_{#tau#tau} (GeV)'
    y_title = 'S/(S+B) weighted events / GeV'
    extra_pad=0.3
    axish = createAxisHists(2,total_bkg,total_bkg.GetXaxis().GetXmin(),total_bkg.GetXaxis().GetXmax()-0.0001)
    axish[1].GetXaxis().SetTitle(x_title)
    axish[1].GetXaxis().SetLabelSize(0.03)
    axish[1].GetXaxis().SetTitleSize(0.04)
    axish[1].GetYaxis().SetNdivisions(4)
    axish[1].GetYaxis().SetTitle("Obs.-Bkg.")
    axish[1].GetYaxis().SetTitleOffset(1.6)
    axish[1].GetYaxis().SetTitleSize(0.04)
    axish[1].GetYaxis().SetLabelSize(0.03)
    axish[0].GetXaxis().SetTitleSize(0)
    axish[0].GetXaxis().SetLabelSize(0)
    axish[0].GetYaxis().SetTitle(y_title)
    axish[0].GetYaxis().SetTitleOffset(1.6)
    axish[0].GetYaxis().SetTitleSize(0.04)
    axish[0].GetYaxis().SetLabelSize(0.03)
    axish[0].SetMinimum(0)
    axish[0].SetMaximum(1.1*(1+extra_pad)*total_bkg.GetMaximum())
    axish[0].Draw()

    stack.Draw("histsame")

    #Draw uncertainty band
    error_hist = total_bkg.Clone()
    error_hist.SetFillColor(CreateTransparentColor(12,0.4))
    error_hist.SetLineColor(CreateTransparentColor(12,0.4))
    error_hist.SetMarkerSize(0)
    error_hist.SetMarkerColor(CreateTransparentColor(12,0.4))
    error_hist.Draw("e2same")

    # draw data
    total_datahist.SetMarkerStyle(20)
    total_datahist.SetLineColor(1)
    #total_datahist.Draw("E same")
    total_datahist.Draw("Esame")
    total_datahist.SetFillStyle(1)
    
    axish[0].Draw("axissame")

    #Setup legend
    legend = plot.PositionedLegend(0.35,0.30,3,0.03)
    legend.SetTextSize(0.025)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)

    legend.AddEntry(total_datahist,"Observed","PEl")
    #Drawn on legend in reverse order looks better
    bkg_histos.reverse()

    background_schemes.reverse()
    leg_hists = [None]*len(bkg_histos)
    for legi,hists in enumerate(bkg_histos):
        legend.AddEntry(hists,background_schemes[legi]['leg_text'],"f")
    total_bkg.SetLineWidth(0)
    legend.AddEntry(error_hist,"Bkg. Uncertainty","f")
    legend.AddEntry(total_sig,"ggH(100 GeV) @ 5.4 pb"%vars(),"l")
    legend.Draw("same")

    # add pT bin label
    bin_label = ''
    if cats[0] == 132: bin_label='p_{T}^{#tau#tau} < 50 GeV'
    if cats[0] == 232: bin_label='50 #leq p_{T}^{#tau#tau} < 100 GeV'
    if cats[0] == 332: bin_label='100 #leq p_{T}^{#tau#tau} < 200 GeV'
    if cats[0] == 432: bin_label='p_{T}^{#tau#tau} #geq 200 GeV'
    if len(cats)==4: bin_label ='p_{T}^{#tau#tau} inclusive'

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)
    latex.DrawLatex(0.21,0.89,"{}".format(bin_label))

    # draw cms and lumi lable
    DrawCMSLogo(pads[0], 'CMS', 'Preliminary', 0, 0.07, -0.0, 2.0, '', 0.85, relExtraDX=0.05)
    lumi = "138 fb^{-1} (13 TeV)"
    plot.DrawTitle(pads[0], lumi, 3, textSize=0.6)


    # draw lower plot with data-bkg

    bkg_sub = SubtractData(error_hist, error_hist)
    data_sub = SubtractData(error_hist, total_datahist)
    if infile2:
      bkgonly.Scale(1.0,"width") 
      bkgonly_sub = SubtractData(error_hist, bkgonly)
    if prefit_plot:
      prefit.Scale(1.0,"width")
      prefit_sub = SubtractData(error_hist, prefit)
    pads[1].cd()
    pads[1].SetGrid(0,1)

    maxy = 0
    miny = 100000
    for i in range(1,data_sub.GetNbinsX()+1):
      hi = data_sub.GetBinContent(i)+data_sub.GetBinError(i)*1.2
      lo = data_sub.GetBinContent(i)-data_sub.GetBinError(i)*1.2
      maxy = max(maxy,hi)
      miny = min(miny,lo)
    for i in range(1,bkg_sub.GetNbinsX()+1):
      hi = bkg_sub.GetBinContent(i)+bkg_sub.GetBinError(i)*1.2
      lo = bkg_sub.GetBinContent(i)-bkg_sub.GetBinError(i)*1.2
      maxy = max(maxy,hi)
      miny = min(miny,lo)
    axish[1].SetMinimum(miny)
    axish[1].SetMaximum(maxy)
    axish[1].Draw("axis")

    if infile2:
      bkgonly_sub.SetLineWidth(2)
      bkgonly_sub.SetLineColor(ROOT.kBlue)
      bkgonly_sub.SetMarkerSize(0)
      bkgonly_sub.SetLineStyle(2)
      bkgonly_sub.Draw("histsame")

      #Setup legend
      legend2 = plot.PositionedLegend(0.25,0.10,3,0.01)
      legend2.SetTextSize(0.025)
      legend2.SetTextFont(42)
      legend2.SetFillStyle(0)

      legend2.AddEntry(bkgonly_sub,"Bkg. only fit","l")
      legend2.AddEntry(total_sig,"Sig.+Bkg. fit","l")
      legend2.Draw("same")

    if prefit_plot:
      prefit_sub.SetLineWidth(3)
      prefit_sub.SetLineColor(ROOT.kMagenta)
      prefit_sub.SetMarkerSize(0)
      prefit_sub.SetLineStyle(2)
      prefit_sub.Draw("histsame")

    total_sig.Draw("histsame")
    data_sub.DrawCopy("e0same")
    bkg_sub.SetMarkerSize(0)
    bkg_sub.SetLineColor(ROOT.kBlack)
    bkg_sub.Draw("e2samehist")

    pads[1].RedrawAxis()


    pads[0].cd()
    #pads[0].GetFrame().Draw()
    pads[0].RedrawAxis()

    plot_name='plot'
    if len(cats) == 4: plot_name+='_cmb'
    elif cats[0] == 432: plot_name+='_GT200'
    elif cats[0] == 332: plot_name+='_100to200'
    elif cats[0] == 232: plot_name+='_50to100'
    elif cats[0] == 132: plot_name+='_0to50'
    c1.SaveAs('prop_plots/%(plot_name)s.pdf' % vars())

HTTMassPlot(cats=[132], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_0to50.root'),infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_0to50.root'))
HTTMassPlot(cats=[232], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_50to100.root'), infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_50to100.root'))
HTTMassPlot(cats=[332], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_100to200.root'), infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_100to200.root'))
HTTMassPlot(cats=[432], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_GT200.root'), infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_GT200.root'))
#HTTMassPlot(cats=[132,232,332,432], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_cmb.root'),infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_cmb.root'))

#HTTMassPlot(cats=[132], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_0to50.root'),infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_0to50.root'),prefit_plot=ROOT.TFile('shapes_prop_plot_prefit_plot_0to50.root'))
#HTTMassPlot(cats=[232], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_50to100.root'), infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_50to100.root'),prefit_plot=ROOT.TFile('shapes_prop_plot_prefit_plot_50to100.root'))
#HTTMassPlot(cats=[332], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_100to200.root'), infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_100to200.root'),prefit_plot=ROOT.TFile('shapes_prop_plot_prefit_plot_100to200.root'))
#HTTMassPlot(cats=[432], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_GT200.root'), infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_GT200.root'),prefit_plot=ROOT.TFile('shapes_prop_plot_prefit_plot_GT200.root'))
#HTTMassPlot(cats=[132,232,332,432], infile=ROOT.TFile('shapes_prop_plot_postfit_plot_cmb.root'),infile2=ROOT.TFile('shapes_prop_plot_postfit_bkgonly_plot_cmb.root'),prefit_plot=ROOT.TFile('shapes_prop_plot_prefit_plot_cmb.root'))
