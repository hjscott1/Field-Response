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

void MakePlot(int planenum, int angle_low, int angle_high, char* runtype1, int binnum1, char* filepath1, char* label1, char* runtype2, int binnum2, char* filepath2, char* label2, char* savefiletext);

void MakePlotsWaveform()
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
      
	  MakePlot(i,angles[k],angles[k+1],(char*)"XX",0,(char*)"./",(char*)"MC",(char*)"XX",0,(char*)"./",(char*)"Data",(char*)"MC_VS_Data");
    }
  }
  
  return;
}

void MakePlot(int planenum, int angle_low, int angle_high, char* runtype1, int binnum1, char* filepath1, char* label1, char* runtype2, int binnum2, char* filepath2, char* label2, char* savefiletext)
{
  TH1F *WaveformHist1;
  TH1F *WaveformHist2;
  
  string bintext1 = "";
  if(binnum1 != 0)
  {
    bintext1 = Form("_bin%d",binnum1);
  }

  string bintext2 = "";
  if(binnum2 != 0)
  {
    bintext2 = Form("_bin%d",binnum2);
  }

  TFile* inputfile_1 = new TFile(Form("./Results/MC/AviWorkflow/100000_Events/WFresults_Plane%d_%s_MC.root",planenum,runtype1),"READ");

  if(angle_high >= angle_low+2)
  {
    WaveformHist1 = (TH1F*) inputfile_1->Get(Form("AnodeRecoHist1D_%dto%d%s",angle_low,angle_low+2,bintext1.c_str()));
  }
  else
  {
    return;
  }

  int angle_1 = angle_low+2;
  while(angle_1 <= angle_high-2)
  {
    WaveformHist1->Add((TH1F*) inputfile_1->Get(Form("AnodeRecoHist1D_%dto%d%s",angle_1,angle_1+2,bintext1.c_str())));
    angle_1 += 2;
  }  
  
  TFile* inputfile_2 = new TFile(Form("./Results/Data/run14608/WFresults_Plane%d_%s_data.root",planenum,runtype2),"READ");
  
  if(angle_high >= angle_low+2)
  {
    WaveformHist2 = (TH1F*) inputfile_2->Get(Form("AnodeRecoHist1D_%dto%d%s",angle_low,angle_low+2,bintext2.c_str()));
  }
  else
  {
    return;
  }

  int angle_2 = angle_low+2;
  while(angle_2 <= angle_high-2)
  {
    WaveformHist2->Add((TH1F*) inputfile_2->Get(Form("AnodeRecoHist1D_%dto%d%s",angle_2,angle_2+2,bintext2.c_str())));
    angle_2 += 2;
  }

  TF1 *smearFunc = new TF1("smearFunc","gaus",timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
  smearFunc->SetParameter(0,1.0);
  smearFunc->SetParameter(1,0.0);
  smearFunc->SetParameter(2,smearFactor);

  TH1F *SmearSimHist = new TH1F("SmearSimHist","",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));

  int numBins = 2*timePadding+1;
  TVirtualFFT::SetTransform(0);

  auto rho_WaveformHist2 = new double[numBins];
  auto ph_WaveformHist2 = new double[numBins];
  TH1 *hm_WaveformHist2 = 0;
  hm_WaveformHist2 = WaveformHist2->FFT(hm_WaveformHist2,"MAG");
  TH1 *hp_WaveformHist2 = 0;
  hp_WaveformHist2 = WaveformHist2->FFT(hp_WaveformHist2,"PH");
  for(int j = 0; j < numBins; j++)
    {
      rho_WaveformHist2[j] = hm_WaveformHist2->GetBinContent(j+1);
      ph_WaveformHist2[j] = hp_WaveformHist2->GetBinContent(j+1);
    }
  delete hm_WaveformHist2;
  delete hp_WaveformHist2;

  for(int i = 1; i <= SmearSimHist->GetNbinsX(); i++)
  {
    SmearSimHist->SetBinContent(i,smearFunc->Eval(SmearSimHist->GetBinCenter(i)));
  }  

  auto rho_SmearSimHist = new double[numBins];
  auto ph_SmearSimHist = new double[numBins];
  TH1 *hm_SmearSimHist = 0;
  hm_SmearSimHist = SmearSimHist->FFT(hm_SmearSimHist,"MAG");
  TH1 *hp_SmearSimHist = 0;
  hp_SmearSimHist = SmearSimHist->FFT(hp_SmearSimHist,"PH");
  for(int j = 0; j < numBins; j++)
  {
    rho_SmearSimHist[j] = hm_SmearSimHist->GetBinContent(j+1);
    ph_SmearSimHist[j] = hp_SmearSimHist->GetBinContent(j+1);
  }
  delete hm_SmearSimHist;
  delete hp_SmearSimHist;

  auto result_re = new double[numBins];
  auto result_im = new double[numBins];

  for(int j = 0; j < numBins; j++)
  {
    float combined_re = (rho_SmearSimHist[j]*cos(ph_SmearSimHist[j])*rho_WaveformHist2[j]*cos(ph_WaveformHist2[j]) - rho_SmearSimHist[j]*sin(ph_SmearSimHist[j])*rho_WaveformHist2[j]*sin(ph_WaveformHist2[j]))/numBins;
    float combined_im = (rho_SmearSimHist[j]*sin(ph_SmearSimHist[j])*rho_WaveformHist2[j]*cos(ph_WaveformHist2[j]) + rho_SmearSimHist[j]*cos(ph_SmearSimHist[j])*rho_WaveformHist2[j]*sin(ph_WaveformHist2[j]))/numBins;

    result_re[j] = combined_re;
    result_im[j] = combined_im;
  }

  TVirtualFFT *ifft;
  ifft = TVirtualFFT::FFT(1,&numBins,"C2R M K");
  ifft->SetPointsComplex(result_re,result_im);
  ifft->Transform();
  
  TH1 *fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");

  TH1F *ResultHist2 = new TH1F("ResultHist2","",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
  
  for(int j = 0; j < numBins; j++)
  {
    ResultHist2->SetBinContent(j+1,fb->GetBinContent(((numBins+j+timePadding) % numBins)+1));
  }

  ResultHist2->Scale(WaveformHist2->Integral()/ResultHist2->Integral());
  
  delete fb;
  delete ifft;
  delete[] rho_WaveformHist2;
  delete[] ph_WaveformHist2;
  delete[] rho_SmearSimHist;
  delete[] ph_SmearSimHist;
  delete[] result_re;
  delete[] result_im;
    
  if((planenum == 0) || (planenum == 1))
  {
    int nbins = ResultHist2->GetNbinsX();
    int minbin_2 = ResultHist2->GetMinimumBin();
    int minbin_1 = WaveformHist1->GetMinimumBin();

    float SF_2 = ResultHist2->GetBinContent(ResultHist2->GetMinimumBin());
    float SF_1 = WaveformHist1->GetBinContent(WaveformHist1->GetMinimumBin());

    WaveformHist1->Scale(SF_2/SF_1);
  }
  else if(planenum == 2)
  {
    WaveformHist1->Scale(ResultHist2->Integral()/WaveformHist1->Integral());
  }
  
  WaveformHist1->SetLineWidth(3.0);
  if(planenum == 0)
  {
    WaveformHist1->SetLineColor(kBlue);
    WaveformHist1->SetMarkerColor(kBlue);
  }
  else if(planenum == 1)
  {
    WaveformHist1->SetLineColor(kGreen+2);
    WaveformHist1->SetMarkerColor(kGreen+2);
  }
  else if(planenum == 2)
  {
    WaveformHist1->SetLineColor(kRed);
    WaveformHist1->SetMarkerColor(kRed);
  }
	    
  ResultHist2->SetLineWidth(3.0);
  ResultHist2->SetLineColor(kBlack);
  ResultHist2->SetMarkerColor(kBlack);

  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.055,"T");
    
  TCanvas c1;
  c1.cd();
  WaveformHist1->Draw("HIST");
  ResultHist2->Draw("HISTsame");
  WaveformHist1->SetTitle(Form("Plane %d, %s TPC, %d\260 - %d\260 tracks: Average Waveform @ Anode",planenum,runtype1,angle_low,angle_high));
  WaveformHist1->GetXaxis()->SetTitle("Time Offset [#mus]");
  WaveformHist1->GetXaxis()->SetTitleSize(0.045);
  WaveformHist1->GetXaxis()->SetTitleOffset(1.05);
  WaveformHist1->GetXaxis()->SetLabelSize(0.04);
  WaveformHist1->GetXaxis()->SetRangeUser(timeTickSF*(-1*max_ticks_plot-0.5),timeTickSF*(max_ticks_plot+0.5));
  //WaveformHist1->GetXaxis()->SetRangeUser(-20, 20);

  if(planenum == 0)
  {
    WaveformHist1->GetYaxis()->SetRangeUser(1.15*min(WaveformHist1->GetMinimum(),ResultHist2->GetMinimum()),max(-0.5*min(WaveformHist1->GetMinimum(),ResultHist2->GetMinimum()),1.2*max(WaveformHist1->GetMaximum(),ResultHist2->GetMaximum())));
  }
  else if(planenum == 1)
  {
    WaveformHist1->GetYaxis()->SetRangeUser(1.2*min(WaveformHist1->GetMinimum(),ResultHist2->GetMinimum()),1.2*max(WaveformHist1->GetMaximum(),ResultHist2->GetMaximum()));
  }
  else if(planenum == 2)
  {
    WaveformHist1->GetYaxis()->SetRangeUser(-0.1*max(WaveformHist1->GetMaximum(),ResultHist2->GetMaximum()),1.1*max(WaveformHist1->GetMaximum(),ResultHist2->GetMaximum()));
  }
  WaveformHist1->GetYaxis()->SetTitle("Arb. Units");
  WaveformHist1->GetYaxis()->SetTitleSize(0.045);
  WaveformHist1->GetYaxis()->SetTitleOffset(1.12);
  WaveformHist1->GetYaxis()->SetLabelSize(0.04);
  TLegend *leg = new TLegend(0.60,0.72,0.85,0.85);
  leg->SetLineColor(kWhite);
  leg->SetTextSize(0.043);
  leg->AddEntry(ResultHist2,label2,"L");
  leg->AddEntry(WaveformHist1,label1,"L");
  leg->Draw("SAME");

  c1.SaveAs(Form("Waveform_Anode_%s_Plane%d_%s_%dto%d.png",savefiletext,planenum,runtype1,angle_low,angle_high));
}
