#include <iostream>
#include <iomanip>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TCut.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVirtualFFT.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"

using namespace std;

//const vector<int> angles = {20, 30, 40, 50, 60, 70, 80};
const vector<int> angles = {20, 80};

const float smearFactor = 0.000001;
//const float smearFactor = 1.0;

//const int max_ticks_plot = 100;
//const int max_ticks_plot = 50;
const int max_ticks_plot = 200;

//const int timePadding = 50;
const int timePadding = 200;
const float timeTickSF = 0.5;

void MakePlot(int planenum, int angle_low, int angle_high, char* runtype, int binnum, char* filepath, char* label);

void MakePlotWaveform()
{
  for(int k = 0; k < angles.size()-1; k++)
  { 
    for(int i = 0; i <= 2; i++)
    {
      char *baseline;
      if(i == 0)
      {
        baseline = (char*) "Right";
      }
      else
      {
        baseline = (char*) "Left";
      }
      
	  MakePlot(i,angles[k],angles[k+1],(char*)"X",0,(char*)"./",(char*)"run14860");
    }
  }
  
  return;
}

void MakePlot(int planenum, int angle_low, int angle_high, char* runtype, int binnum, char* filepath, char* label)
{
  TH1F *WaveformHist;
  
  string bintext = "";
  if(binnum != 0)
  {
    bintext = Form("_bin%d",binnum);
  }

  TFile* inputfile = new TFile(Form("./Results/data/run14860/WFresults_Plane%d_%s_data.root",planenum,runtype),"READ");

  if(angle_high >= angle_low+2)
  {
    WaveformHist = (TH1F*) inputfile->Get(Form("AnodeRecoHist1D_%dto%d%s",angle_low,angle_low+2,bintext.c_str()));
  }
  else
  {
    return;
  }

  int angle = angle_low+2;
  while(angle <= angle_high-2)
  {
    WaveformHist->Add((TH1F*) inputfile->Get(Form("AnodeRecoHist1D_%dto%d%s",angle,angle+2,bintext.c_str())));
    angle += 2;
  }  
  
  WaveformHist->SetLineWidth(3.0);
  if(planenum == 0)
  {
    WaveformHist->SetLineColor(kBlue);
    WaveformHist->SetMarkerColor(kBlue);
  }
  else if(planenum == 1)
  {
    WaveformHist->SetLineColor(kGreen+2);
    WaveformHist->SetMarkerColor(kGreen+2);
  }
  else if(planenum == 2)
  {
    WaveformHist->SetLineColor(kRed);
    WaveformHist->SetMarkerColor(kRed);
  }
	    
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.055,"T");
    
  TCanvas c1;
  c1.cd();
  WaveformHist->Draw("HIST");
  WaveformHist->SetTitle(Form("Plane %d, %d\260 - %d\260 tracks: Average Waveform @ Anode",planenum,angle_low,angle_high));
  WaveformHist->GetXaxis()->SetTitle("Time Offset [#mus]");
  WaveformHist->GetXaxis()->SetTitleSize(0.045);
  WaveformHist->GetXaxis()->SetTitleOffset(1.05);
  WaveformHist->GetXaxis()->SetLabelSize(0.04);
  WaveformHist->GetXaxis()->SetRangeUser(timeTickSF*(-1*max_ticks_plot-0.5),timeTickSF*(max_ticks_plot+0.5));
  //WaveformHist1->GetXaxis()->SetRangeUser(-20, 20);

  WaveformHist->GetYaxis()->SetTitle("Arb. Units");
  WaveformHist->GetYaxis()->SetTitleSize(0.045);
  WaveformHist->GetYaxis()->SetTitleOffset(1.12);
  WaveformHist->GetYaxis()->SetLabelSize(0.04);
  TLegend *leg = new TLegend(0.60,0.72,0.85,0.85);
  leg->SetLineColor(kWhite);
  leg->SetTextSize(0.043);
  leg->AddEntry(WaveformHist,label,"L");
  leg->Draw("SAME");

  c1.SaveAs(Form("Results/data/run14860/Waveform_Anode_Plane%d_%s_%dto%d.png",planenum,runtype,angle_low,angle_high));
}
