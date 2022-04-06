#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import CombineHarvester.CombineTools.plotting as plot
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    'input', nargs='+', help="""Input json files""")
parser.add_argument(
    '--output', '-o', default='limit', help="""Name of the output
    plot without file extension""")
parser.add_argument(
    '--show', default='exp,obs')
parser.add_argument(
    '--x-title', default='m_{#phi} (GeV)', help="""Title for the x-axis""")
parser.add_argument(
    '--y-title', default=None, help="""Title for the y-axis""")
parser.add_argument(
    '--y-axis-min', default=None, help="""Minimum for y-axis range""")
parser.add_argument(
    '--y-axis-max', default=None, help="""Maximum for y-axis range""")
parser.add_argument(
    '--process', choices=['gg#phi','bb#phi'], help='The process on which a limit has been calculated.', default="gg#phi")
parser.add_argument(
    '--cms-sub', default='Internal', help="""Text below the CMS logo""")
parser.add_argument(
    '--title-right', default='', help="""Right header text above the frame""")
parser.add_argument(
    '--title-left', default='', help="""Left header text above the frame""")
parser.add_argument(
    '--logy', action='store_true', help="""Draw y-axis in log scale""")
parser.add_argument(
    '--logx', action='store_true', help="""Draw x-axis in log scale""")
parser.add_argument(
    '--low_high_split', action='store_true', help="""Draw line and labels for low mass - high mass split""")
parser.add_argument(
    '--ratio-to', default=None)
parser.add_argument(
    '--pad-style', default=None, help="""Extra style options for the pad, e.g. Grid=(1,1)""")
parser.add_argument(
    '--auto-style', nargs='?', const='', default=None, help="""Take line colors and styles from a pre-defined list""")
args = parser.parse_args()

style_dict_hig_17_020 = {
        'style' : {
            'exp0' : { 'LineColor' : ROOT.kBlack, 'LineStyle' : 2},
            'exp1' : { 'FillColor' : ROOT.kGreen+1}, 
            'exp2' : { 'FillColor' : ROOT.kOrange}
            },
        'legend' : {
            'exp1' : { 'Label' : '68% expected'},
            'exp2' : { 'Label' : '95% expected'}
            }

        }

style_dict = style_dict_hig_17_020

def DrawAxisHists(pads, axis_hists, def_pad=None):
    for i, pad in enumerate(pads):
        pad.cd()
        axis_hists[i].Draw('AXIS')
        axis_hists[i].Draw('AXIGSAME')
    if def_pad is not None:
        def_pad.cd()

## Boilerplate
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()
ROOT.gStyle.SetNdivisions(510, 'XYZ') # probably looks better

canv = ROOT.TCanvas(args.output, args.output)

if args.ratio_to is not None:
    pads = plot.TwoPadSplit(0.30, 0.01, 0.01)
else:
    pads = plot.OnePad()

# Set the style options of the pads
for padx in pads:
    # Use tick marks on oppsite axis edges
    plot.Set(padx, Tickx=1, Ticky=1, Logx=args.logx)
    if args.pad_style is not None:
        settings = {x.split('=')[0]: eval(x.split('=')[1]) for x in args.pad_style.split(',')}
        print 'Applying style options to the TPad(s):'
        print settings
        plot.Set(padx, **settings)

graphs = []
graph_sets = []

legend = plot.PositionedLegend(0.48, 0.25, 3, 0.015)
legend.SetTextSize(0.04)

axis = None

defcols = [
    ROOT.kGreen+3, ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kYellow+2,
    ROOT.kOrange+10, ROOT.kCyan+3, ROOT.kMagenta+2, ROOT.kViolet-5, ROOT.kGray
    ]

deflines = [1, 2, 3]

if args.auto_style is not None:
    icol = {x: 0 for x in args.auto_style.split(',')}
    icol['default'] = 0
    iline = {}
    iline['default'] = 1
    for i, x in enumerate(args.auto_style.split(',')):
        iline[x] = i+1

# Process each input argument
has_band = False

dummyhist = ROOT.TH1F("dummy", "", 1, 0, 1)
plot.Set(dummyhist, LineColor=ROOT.kWhite, FillColor=ROOT.kWhite)

def RemovePoints(graph_set, high=True):
  graph_set_new = {}
  for key,g in graph_set.items(): 
    x=ROOT.Double()
    y=ROOT.Double()
    
    import math
    shift=5.
    g_clone = g.Clone()
    for i in range(g_clone.GetN(),-1,-1):
      g_clone.GetPoint(i,x,y)
      if high and key =='obs' and x < 250: g_clone.RemovePoint(i)
      elif high and key !='obs' and x < 250:
        x_bound=225.
        x_bound+=1.0*shift
        #x_bound=200
        if x==200:
          y_new = g_clone.Eval(x_bound)
          if key in ['exp1','exp2']:
            print key
            x_=ROOT.Double()
            y_=ROOT.Double()
            g.GetPoint(i+1,x_,y_)

            lo_0 = y - g.GetErrorYlow(i)
            lo_1 = y_ - g.GetErrorYlow(i+1)
            hi_0 = y + g.GetErrorYhigh(i)
            hi_1 = y_ + g.GetErrorYhigh(i+1)

            up_new=(hi_1-hi_0)/(x_-x)*(x_bound-x) + hi_0 -y_new
            down_new = y_new - ((lo_1-lo_0)/(x_-x)*(x_bound-x)+lo_0)

#            down_new = y_new - ((lo_1-lo_0)/(x-x_)*(x_bound-x_)+lo_0)
            #up_new = ((hi_1-hi_0)/(x-x_)*(x_bound-x_)+hi_0) -y_new

            g_clone.SetPointEYhigh(i,up_new)
            g_clone.SetPointEYlow(i,down_new)
          g_clone.SetPoint(i,x_bound,y_new) 
        else: g_clone.RemovePoint(i)
      if not high and key == 'obs' and x>= 250: g_clone.RemovePoint(i)
      elif not high and key != 'obs' and x> 200: 
        x_bound=225.-shift
        if x==250:
          y_new = g.Eval(x_bound)
          if key in ['exp1','exp2']:
            x0=ROOT.Double()
            y0=ROOT.Double()
            g.GetPoint(i-1,x0,y0)

            lo_0 = y0 - g.GetErrorYlow(i-1) 
            lo_1 = y - g.GetErrorYlow(i)
            hi_0 = y0 + g.GetErrorYhigh(i-1)
            hi_1 = y + g.GetErrorYhigh(i)
            down_new = y_new - ((lo_1-lo_0)/(x-x0)*(x_bound-x0)+lo_0)
            up_new = ((hi_1-hi_0)/(x-x0)*(x_bound-x0)+hi_0) -y_new

            g_clone.SetPointEYhigh(i,up_new)
            g_clone.SetPointEYlow(i,down_new)
          g_clone.SetPoint(i,x_bound,y_new)
        else: g_clone.RemovePoint(i) 


    graph_set_new[key] = g_clone
  return graph_set_new

def DrawLimitBandWithRange(pad, graph_dict, draw=['exp2', 'exp1', 'exp0', 'obs'], draw_legend=None,
                  legend=None, legend_overwrite=None, range_=None):
    legend_dict = {
        'obs' : { 'Label' : 'Observed', 'LegendStyle' : 'LP', 'DrawStyle' : 'PLSAME'},
        'exp0' : { 'Label' : 'Expected', 'LegendStyle' : 'L', 'DrawStyle' : 'LSAME'},
        'exp1' : { 'Label' : '#pm1#sigma Expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'},
        'exp2' : { 'Label' : '#pm2#sigma Expected', 'LegendStyle' : 'F', 'DrawStyle' : '3SAME'}
    }
    if legend_overwrite is not None:
        for key in legend_overwrite:
            if key in legend_dict:
                legend_dict[key].update(legend_overwrite[key])
            else:
                legend_dict[key] = legend_overwrite[key]
    pad.cd()
    for key in draw:
        if key in graph_dict:
            graph_dict[key].Draw(legend_dict[key]['DrawStyle'])
    if legend is not None:
        if draw_legend is None:
            draw_legend = reversed(draw)
        for key in draw_legend:
            if key in graph_dict:
                legend.AddEntry(graph_dict[key],legend_dict[key]['Label'],legend_dict[key]['LegendStyle'])


for i, src in enumerate(args.input):
    splitsrc = src.split(':')
    file = splitsrc[0]
    # limit.json => Draw as full obs + exp limit band
    if len(splitsrc) == 1:
        graph_sets.append(plot.StandardLimitsFromJSONFile(file, args.show.split(',')))
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), graph_sets[-1].values()[0], True)
            for a in axis: a.GetXaxis().SetLimits(60., 3500,)
            DrawAxisHists(pads, axis, pads[0])
        plot.StyleLimitBand(graph_sets[-1],overwrite_style_dict=style_dict["style"])
         
        if not args.low_high_split: plot.DrawLimitBand(pads[0], graph_sets[-1], legend=legend,legend_overwrite=style_dict["legend"])
        else: 
          if i==1: 
            graph_set_high = RemovePoints(graph_sets[-1], high=True)
            DrawLimitBandWithRange(pads[0], graph_set_high, legend=legend,legend_overwrite=style_dict["legend"],range_=None)
          if i==0: 
            graph_set_low = RemovePoints(graph_sets[-1], high=False)
            DrawLimitBandWithRange(pads[0], graph_set_low, range_=None)
        pads[0].RedrawAxis()
        pads[0].RedrawAxis('g')
        pads[0].GetFrame().Draw()
        has_band = True  # useful to know later if we want to do style settings
                         # based on whether or not the expected band has been drawn

        #break

    # limit.json:X => Draw a single graph for entry X in the json file 
    # 'limit.json:X:Title="Blah",LineColor=4,...' =>
    # as before but also apply style options to TGraph

    elif len(splitsrc) >= 2:
        settings = {}
        settings['Title'] = src
        if args.auto_style is not None:
            nm = 'default'
            for x in icol.keys():
                if x in splitsrc[1]:
                    nm = x
            i = icol[nm]  # take the next default color...
            j = iline[nm]  # take the next default line style...
            settings['LineColor'] = defcols[i]
            settings['MarkerColor'] = defcols[i]
            settings['LineStyle'] = j
            icol[nm] = (i+1) if (i+1) < len(defcols) else 0
        graph = plot.LimitTGraphFromJSONFile(file, splitsrc[1])
        graphs.append(graph)
        if len(splitsrc) >= 3:
            settings.update({x.split('=')[0]: eval(x.split('=')[1]) for x in splitsrc[2].split(',')})
        plot.Set(graphs[-1], **settings)
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), graphs[-1], True)
            DrawAxisHists(pads, axis, pads[0])
        graphs[-1].Draw('PLSAME')
        legend.AddEntry(graphs[-1], '', 'PL')



axis[0].GetYaxis().SetTitle('95% CL limit on #sigma#font[42]{(gg#phi)}#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)')
if args.process == "bb#phi":
    axis[0].GetYaxis().SetTitle('95% CL limit on #sigma#font[42]{(bb#phi)}#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)')
if args.y_title is not None:
    axis[0].GetYaxis().SetTitle(args.y_title)
axis[0].GetXaxis().SetTitle(args.x_title)
axis[0].GetXaxis().SetNoExponent()
axis[0].GetXaxis().SetMoreLogLabels()
axis[0].GetXaxis().SetLabelOffset(axis[0].GetXaxis().GetLabelOffset()*2)

if args.logy:
    axis[0].SetMinimum(0.1)  # we'll fix this later
    pads[0].SetLogy(True)
    # Apparently switching to logy puts the band back over the top of the axis
    if has_band:
        pads[0].RedrawAxis()
        pads[0].RedrawAxis('g')
        pads[0].GetFrame().Draw()
    # axis[0].GetYaxis().SetMoreLogLabels()
    # axis[0].SetNdivisions(50005, "X")

y_min, y_max = (plot.GetPadYMin(pads[0]), plot.GetPadYMax(pads[0]))
plot.FixBothRanges(pads[0], y_min if args.logy else 0, 0.05 if args.logy else 0, y_max, 0.25)

if args.y_axis_min is not None or args.y_axis_max is not None:
  hobj = plot.GetAxisHist(pads[0])
  if args.y_axis_min is not None: hobj.SetMinimum(float(args.y_axis_min))
  if args.y_axis_max is not None: hobj.SetMaximum(float(args.y_axis_max))


ratio_graph_sets = []
ratio_graphs = []

if args.ratio_to is not None:
    pads[1].cd()
    plot.SetupTwoPadSplitAsRatio(pads, axis[0], axis[1], '', True, 0.1, 2.4)
    axis[1].SetNdivisions(506, 'Y')
    splitsrc = args.ratio_to.split(':')
    ref = plot.LimitTGraphFromJSONFile(splitsrc[0], splitsrc[1])
    for gr_set in graph_sets:
        ratio_set = {}
        for key in gr_set:
            ratio_set[key] = plot.GraphDivide(gr_set[key], ref)
        ratio_graph_sets.append(ratio_set)
        plot.DrawLimitBand(pads[1], ratio_graph_sets[-1])
        pads[1].RedrawAxis()
        pads[1].RedrawAxis('g')
        pads[1].GetFrame().Draw()
    for gr in graphs:
        ratio_graphs.append(plot.GraphDivide(gr, ref))
        ratio_graphs[-1].Draw('LP')
    ry_min, ry_max = (plot.GetPadYMin(pads[1]), plot.GetPadYMax(pads[1]))
    plot.FixBothRanges(pads[1], ry_min, 0.1, ry_max, 0.1)


pads[0].cd()
if legend.GetNRows() == 1:
    legend.SetY1(legend.GetY2() - 0.5*(legend.GetY2()-legend.GetY1()))
legend.Draw()

plot.DrawTitle(pads[0], args.title_right % vars(), 3)
#plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 1, 0.045, 0.05, 1.0, '', 0.9)
plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 0, 0.15, 0, 0, '', 0.9)
#  plot.DrawCMSLogo(c, 'CMS', 'Preliminary', 0, 0.15, 0, 0, '', 0.9)

#plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 11, 0.045, 0.035, 1.2, '', 0.8)
#plot.DrawTitle(pads[0], args.title_right, 3)
#plot.DrawTitle(pads[0], args.title_left, 1)


latex2 = ROOT.TLatex()
latex2.SetNDC()
latex2.SetTextAngle(0)
latex2.SetTextAlign(12)
latex2.SetTextFont(42)
latex2.SetTextSize(0.04)

if args.low_high_split:

  line =  ROOT.TLine(225.,hobj.GetMinimum(),225,hobj.GetMaximum())
  line.Draw()
  latex2.DrawLatex(0.19,0.17, 'Low-mass')
  latex2.DrawLatex(0.45,0.17, 'High-mass')
  #latex2.SetTextSize(0.05)

if 'freezebbH' in args.output:
  latex2.DrawLatex(0.48,0.62, 'bb#phi set to zero')

if 'freezeggH' in args.output:
  latex2.DrawLatex(0.48,0.62, 'gg#phi set to zero')

if 'Tonly' in args.output:
  latex2.DrawLatex(0.48,0.62, 't quark only')

if 'Bonly' in args.output:
  latex2.DrawLatex(0.48,0.62, 'b quark only')

canv.Print('.pdf')
canv.Print('.png')
