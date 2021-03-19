#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ROOT
import CombineHarvester.CombineTools.plotting as plot
import argparse

class Color(int):
    """Create a new ROOT.TColor object with an associated index"""
    __slots__ = ["object", "name"]

    def __new__(cls, r, g, b, name=""):
        self = int.__new__(cls, ROOT.TColor.GetFreeColorIndex())
        self.object = ROOT.TColor(self, r, g, b, name, 1.0)
        self.name = name
        return self


parser = argparse.ArgumentParser()
parser.add_argument(
    'input', nargs='+', help="""Input json files""")
parser.add_argument(
    '--output', '-o', default='limit', help="""Name of the output
    plot without file extension""")
parser.add_argument(
    '--show', default='exp,obs')
parser.add_argument(
    '--x-title', default="m_{h_{S}} (GeV)", help="""Title for the x-axis""")
parser.add_argument(
    '--y-title', default=None, help="""Title for the y-axis""")
parser.add_argument(
    '--y-axis-min', default=None, help="""Minimum for y-axis range""")
parser.add_argument(
    '--y-axis-max', default=None, help="""Maximum for y-axis range""")
parser.add_argument(
    '--process', choices=['gg#phi','bb#phi','nmssm'], help='The process on which a limit has been calculated.', default="gg#phi")
parser.add_argument(
    '--cms-sub', default='', help="""Text below the CMS logo""")
parser.add_argument(
    '--title-right', default='', help="""Right header text above the frame""")
parser.add_argument(
    '--title-left', default='', help="""Left header text above the frame""")
parser.add_argument(
    '--logy', action='store_true', help="""Draw y-axis in log scale""")
parser.add_argument(
    '--logx', action='store_true', help="""Draw x-axis in log scale""")
parser.add_argument(
    '--ratio-to', default=None)
parser.add_argument(
    '--xmax', default=None, type=float)
parser.add_argument(
    '--mass', default=None, type=int)
parser.add_argument(
    '--pad-style', default=None, help="""Extra style options for the pad, e.g. Grid=(1,1)""")
parser.add_argument(
    '--auto-style', nargs='?', const='', default=None, help="""Take line colors and styles from a pre-defined list""")
args = parser.parse_args()

# style_dict_hig_17_020 = {
#         'style' : {
#             'exp0' : { 'LineColor' : ROOT.kBlack, 'LineStyle' : 2},
#             'exp1' : { 'FillColor' : ROOT.kGreen+1}, 
#             'exp2' : { 'FillColor' : ROOT.kOrange}
#             },
#         'legend' : {
#             'exp1' : { 'Label' : '68% expected'},
#             'exp2' : { 'Label' : '95% expected'}
#             }

#         }
colors = [
    Color(0.87, 0.73, 0.53, "myyellow"),
    Color(0.57, 0.69, 0.32, "mygreen"),
    Color(0.93, 0.65, 0.17, "myorange"),
    Color(0.96, 0.65, 0.36, "myorange3"),
    Color(0.97, 0.91, 0.40, "myyellow3"),
    Color(0.66, 0.81, 0.33, "mygreen3"),
                ]
                
for color in colors:
    setattr(ROOT, color.name, color)

    
style_dict_hig_17_020 = {
        'style' : {
            'exp0' : { 'LineColor' : 12, 'LineStyle' : 7, 'LineWidth' : 3},
            'exp1' : { 'FillColor' : ROOT.mygreen3}, 
            'exp2' : { 'FillColor' : ROOT.myyellow3}
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
    pads = plot.TwoPadSplit(0.4, 0.01, 0.01)
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

legend = plot.PositionedLegend(0.24, 0.2, 3, 0.015)
legend.SetTextSize(0.03)

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
scale_factors = {
    "240": 22,
    "280": 21,
    "320": 20,
    "360": 19,
    "400": 18,
    "450": 17,
    "500": 16,
    "550": 15,
    "600": 14,
    "700": 13,
    "800": 12,
    "900": 11,
    "1000": 10,
    "1200": 9,
    "1400": 8,
    "1600": 7,
    "1800": 6,
    "2000": 4,
    "2500": 2,
    "3000": 0

}
texts = []
arrows = []
xtemp, ytemp = ROOT.Double(0.0), ROOT.Double(0.0)

for i,src in enumerate(args.input):
    mass =  src.split("_")[4]
    splitsrc = src.split(':')
    file = splitsrc[0]
    # limit.json => Draw as full obs + exp limit band
    if len(splitsrc) == 1:
        graph_sets.append(plot.StandardLimitsFromJSONFile(file, args.show.split(',')))
        if axis is None:
            axis = plot.CreateAxisHists(len(pads), graph_sets[-1].values()[0], True)
            DrawAxisHists(pads, axis, pads[0])
        plot.StyleLimitBand(graph_sets[-1],overwrite_style_dict=style_dict["style"])
        filler = ""
        # style_dict["legend"]["obs"] = {"Label": "m_{H} = %s GeV%s (#times 10^{%s})" % (mass, filler, str(scale_factors[mass]))}
        if i==0:
            plot.DrawLimitBand(pads[0], graph_sets[-1], legend=legend,legend_overwrite=style_dict["legend"])
        else:
            plot.DrawLimitBand(pads[0], graph_sets[-1])
        if True:
            if not mass in ["2000","2500","3000"]:
                texts.append(ROOT.TLatex(1.12*(float(mass)),10**(scale_factors[mass]-1),"m_{H} = %s GeV%s (#times 10^{%s})" % (mass, filler, str(scale_factors[mass]))))
                texts[-1].SetTextFont(42)
                texts[-1].SetTextSize(0.025)
                texts[-1].Draw()

                graph_sets[-1]["exp0"].GetPoint(graph_sets[-1]["exp0"].GetN()-1,xtemp,ytemp)

                arrows.append(ROOT.TArrow(1.12*(float(mass)-10),10**(scale_factors[mass]-1),1.05*xtemp,1.05*ytemp,0.01,"|>"))
                arrows[-1].Draw()
            elif mass=="2000":
                texts.append(ROOT.TLatex(1.1*(float(mass)),10**(scale_factors[mass]-1),"m_{H} = %s GeV%s (#times 10^{%s})" % (mass, filler, str(scale_factors[mass]))))
                texts[-1].SetTextFont(42)
                texts[-1].SetTextSize(0.025)
                texts[-1].Draw()

                graph_sets[-1]["exp0"].GetPoint(graph_sets[-1]["exp0"].GetN()-1,xtemp,ytemp)

                arrows.append(ROOT.TArrow(1.1*(float(mass)),10**(scale_factors[mass]-1),1.05*xtemp,1.05*ytemp,0.01,"|>"))
                arrows[-1].Draw()
            elif mass=="2500":
                texts.append(ROOT.TLatex(1.0*(float(mass)-340),10**(scale_factors[mass]-1),"m_{H} = %s GeV%s (#times 10^{%s})" % (mass, filler, str(scale_factors[mass]))))
                texts[-1].SetTextFont(42)
                texts[-1].SetTextSize(0.024)
                texts[-1].Draw()

                graph_sets[-1]["exp0"].GetPoint(graph_sets[-1]["exp0"].GetN()-1,xtemp,ytemp)

                arrows.append(ROOT.TArrow(1.0*(float(mass)-340),10**(scale_factors[mass]-1),0.9*xtemp,2.0*ytemp,0.01,"|>"))
                arrows[-1].Draw()                
            else:
                texts.append(ROOT.TLatex(int(mass)-600,1.1*10**(scale_factors[mass]-1),"m_{H} = %s GeV%s (#times 10^{%s})" % (mass, filler, str(scale_factors[mass]))))
                texts[-1].SetTextFont(42)
                texts[-1].SetTextSize(0.0233)
                texts[-1].Draw()

                graph_sets[-1]["exp0"].GetPoint(graph_sets[-1]["exp0"].GetN()-1,xtemp,ytemp)

                arrows.append(ROOT.TArrow(int(mass)-600,10**(scale_factors[mass]-1),2300.,0.01,0.01,"|>"))
                arrows[-1].Draw()            
        pads[0].RedrawAxis()
        pads[0].RedrawAxis('g')
        pads[0].GetFrame().Draw()

        has_band = True  # useful to know later if we want to do style settings
                         # based on whether or not the expected band has been drawn

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

mass = args.mass

axis[0].GetYaxis().SetTitle('95% CL limit on #sigma#font[42]{(gg#phi)}#upoint#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)')
if args.process == "bb#phi":
    axis[0].GetYaxis().SetTitle('95% CL limit on #sigma#font[42]{(bb#phi)}#upoint#font[52]{B}#font[42]{(#phi#rightarrow#tau#tau)} (pb)')
if args.process == "nmssm":
    axis[0].GetYaxis().SetTitle("#scale[0.9]{95% CL limit on #sigma#times #font[52]{B}#font[42]{(H#rightarrowh(#tau#tau) h_{S}(bb))}  (pb)}")
if args.y_title is not None:
    axis[0].GetYaxis().SetTitle(args.y_title)
axis[0].GetXaxis().SetTitle(args.x_title)
axis[0].GetXaxis().SetNoExponent()
axis[0].GetXaxis().SetMoreLogLabels()
axis[0].GetXaxis().SetLabelOffset(axis[0].GetXaxis().GetLabelOffset()*2)
if args.xmax is not None:
    axis[0].GetXaxis().SetLimits(59.9,args.xmax)

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

xbatches = [85,110,170,350,550,750] if mass<1001. else [120,350,800,1300,2000]

batch_lines = []
for x_batch in xbatches:
    batch_line = ROOT.TLine(x_batch,float(args.y_axis_min),x_batch,float(args.y_axis_max))
    batch_line.SetLineColor(12)
    batch_line.SetLineStyle(2)
    batch_lines.append(batch_line)
# for batch_line in batch_lines:        
#     batch_line.Draw("SAME")

pads[0].cd()
if legend.GetNRows() == 1:
    legend.SetY1(legend.GetY2() - 0.5*(legend.GetY2()-legend.GetY1()))
legend.SetFillStyle(0)
legend.Draw()

channel_label = {"mt": "#mu^{}_{}#tau^{}_{h}",
                "tt": "#tau^{}_{h}#tau^{}_{h}",
                "et":  "e^{}_{}#tau^{}_{h}",
                "em": "e#mu",
                "all": "e^{}_{}#tau^{}_{h}+#mu^{}_{}#tau^{}_{h}+#tau^{}_{h}#tau^{}_{h}"
                }
for ch in channel_label.keys():
    if ch in args.title_left:
        title_left = args.title_left.replace(ch,channel_label[ch])
        break
plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 11, 0.045, 0.035, 1.2, '', 0.8)
plot.DrawTitle(pads[0], args.title_right, 3)
plot.DrawTitle(pads[0], title_left, 1)

canv.Print('.pdf')
canv.Print('.png')
