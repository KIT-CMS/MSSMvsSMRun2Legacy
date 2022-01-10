#include <TH1.h>
#include <TPad.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TROOT.h>

/// define common styles
void SetStyle();
/// create a canveas in common style (default is square)
TCanvas* MakeCanvas(const char* name, const char* title, int dX=600, int dY=600);
/// define common style of pads
void InitSubPad(TPad* pad, int i);
/// define a histogram file (not yet including poly markers, which are used fro data)
void InitHist(TH1* hist, const char* xtit, const char* ytit="Number of Entries", int color=kBlack, int style=0);
/// define poly markers and line width for data
void InitData(TH1* hist);
/// define line style, width and color for signal sample
void InitSignal(TH1* hist);
/// define common legend style 
void SetLegendStyle(TLegend* leg);
/// add CMSPreliminary logo, lumi and channel
void CMSPrelim(const char* dataset, double lowX=0.4, double lowY=0.74);
void plotchannel(TString channel);



/*
  Implementation
*/


void SetStyle()
{
  TGaxis::SetMaxDigits(4);
  TStyle *HttStyle = new TStyle("Htt-Style","The Perfect Style for Plots ;-)");
  gStyle = HttStyle;

  // Canvas
  HttStyle->SetCanvasColor     (0);
  HttStyle->SetCanvasBorderSize(10);
  HttStyle->SetCanvasBorderMode(0);
  HttStyle->SetCanvasDefH      (700);
  HttStyle->SetCanvasDefW      (700);
  HttStyle->SetCanvasDefX      (100);
  HttStyle->SetCanvasDefY      (100);

  // color palette for 2D temperature plots
  HttStyle->SetPalette(1,0);

  // Pads
  HttStyle->SetPadColor       (0);
  HttStyle->SetPadBorderSize  (10);
  HttStyle->SetPadBorderMode  (0);
  HttStyle->SetPadBottomMargin(0.65);
  HttStyle->SetPadTopMargin   (0.06);
  HttStyle->SetPadLeftMargin  (0.09);
  HttStyle->SetPadRightMargin (0.06);
  HttStyle->SetPadGridX       (0);
  HttStyle->SetPadGridY       (0);
  HttStyle->SetPadTickX       (1);
  HttStyle->SetPadTickY       (1);

  // Frames
  HttStyle->SetLineWidth(3);
  HttStyle->SetFrameFillStyle ( 0);
  HttStyle->SetFrameFillColor ( 0);
  HttStyle->SetFrameLineColor ( 1);
  HttStyle->SetFrameLineStyle ( 0);
  HttStyle->SetFrameLineWidth ( 2);
  HttStyle->SetFrameBorderSize(10);
  HttStyle->SetFrameBorderMode( 0);

  // Histograms
  HttStyle->SetHistFillColor(2);
  //HttStyle->SetHistFillStyle(0);
  HttStyle->SetHistLineColor(1);
  HttStyle->SetHistLineStyle(0);
  HttStyle->SetHistLineWidth(2);
  HttStyle->SetNdivisions(506,"XY");

  // Functions
  HttStyle->SetFuncColor(1);
  HttStyle->SetFuncStyle(0);
  HttStyle->SetFuncWidth(2);

  // Various
  HttStyle->SetMarkerStyle(20);
  HttStyle->SetMarkerColor(kBlack);
  HttStyle->SetMarkerSize (1.4);

  HttStyle->SetTitleBorderSize(0);
  HttStyle->SetTitleFillColor (0);
  //HttStyle->SetTitleX         (0.2);

  HttStyle->SetStripDecimals(kFALSE);
  HttStyle->SetLineStyleString(11,"20 10");

  HttStyle->SetTitleSize  (30,"XY");
  HttStyle->SetTitleOffset(1.60,"X");
  HttStyle->SetTitleOffset(2.10,"Y");
  HttStyle->SetLabelOffset(0.010,"XY");
  HttStyle->SetLabelFont  (43,"XY");
  HttStyle->SetLabelSize  (30,"XY");

  //HttStyle->SetTextFont   (43);
  //HttStyle->SetTextSize   (30);

  HttStyle->SetStatFont   (43);
  HttStyle->SetTitleFont  (43);
  HttStyle->SetTitleFont  (43,"XY");

  HttStyle->SetOptStat    (0);

  gROOT->ForceStyle();
  return;
}

TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.09);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.65);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void InitSubPad(TPad* pad, int i)
{
  pad->cd(i);
  TPad *tmpPad = (TPad*) pad->GetPad(i);
  tmpPad->SetLeftMargin  (0.09);
  tmpPad->SetTopMargin   (0.05);
  tmpPad->SetRightMargin (0.07);
  tmpPad->SetBottomMargin(0.65);
  return;
}

void InitSignal(TH1 *hist)
{
  hist->SetFillStyle(0.);
  hist->SetLineStyle(11);
  hist->SetLineWidth(2.);
  hist->SetLineColor(kBlue+3);
}

void InitHist(TH1 *hist, const char *xtit, const char *ytit, int color, int style)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(    2.);
  hist->SetFillColor(color );
  hist->SetFillStyle(style );
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.200,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.040,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.300,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.050,"X");
  hist->SetLabelFont  (42   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (0.6);
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->SetTitle("");  
  return;
}

void InitData(TH1* hist)
{
  hist->SetMarkerStyle(20.);
  hist->SetMarkerSize (1.3);
  hist->SetLineWidth  ( 2.);
}

void SetLegendStyle(TLegend* leg) 
{
  leg->SetFillStyle (0);
  leg->SetFillColor (0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
}

void CMSPrelim(const char* dataset, double lowX, double lowY)
{
  TPaveText* cmsprel  = new TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC");
  cmsprel->SetBorderSize(   0 );
  cmsprel->SetFillStyle(    0 );
  cmsprel->SetTextAlign(   12 );
  cmsprel->SetTextSize ( 0.05 );
  cmsprel->SetTextColor(    1 );
  cmsprel->SetTextFont (   62 );
  cmsprel->AddText("CMS");
  cmsprel->Draw();

  TPaveText* lumi     = new TPaveText(lowX+0.08, lowY+0.061, lowX+0.45, lowY+0.161, "NDC");
  lumi->SetBorderSize(   0 );
  lumi->SetFillStyle(    0 );
  lumi->SetTextAlign(   12 );
  lumi->SetTextSize ( 0.05 );
  lumi->SetTextColor(    1 );
  lumi->SetTextFont (   62 );
  lumi->AddText(dataset);
  lumi->Draw();

}
void plotchannel(TString channel) {

   TLatex * tex = new TLatex(0.2,0.94,channel);
   tex->SetNDC();
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();

}
