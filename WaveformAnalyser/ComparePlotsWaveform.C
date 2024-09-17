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

//const vector<int> angles = {20, 40, 60, 80};
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

void ComparePlotsWaveform()
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

			MakePlot(i,angles[k],angles[k+1],(char*)"E",0,(char*)"./",(char*)"E",(char*)"W",0,(char*)"./",(char*)"W",(char*)"E_VS_W");
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

	TFile* inputfile_1 = new TFile(Form("./Results/data/WireBias2/WFresults_Plane%d_%s_data.root",planenum,runtype1), "READ");
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
	
	TFile* inputfile_2 = new TFile(Form("./Results/data/WireBias2/WFresults_Plane%d_%s_data.root",planenum,runtype2), "READ");
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

	WaveformHist1->SetLineWidth(3.0);
	WaveformHist2->SetLineWidth(3.0);
	if(planenum == 0)
	{
		WaveformHist1->SetLineColor(kBlue);
		WaveformHist1->SetMarkerColor(kBlue);
		WaveformHist2->SetLineColor(kBlack);
		WaveformHist2->SetMarkerColor(kBlack);
	}
	else if(planenum == 1)
	{
		WaveformHist1->SetLineColor(kGreen+2);
		WaveformHist1->SetMarkerColor(kGreen+2);
		WaveformHist2->SetLineColor(kBlack);
		WaveformHist2->SetMarkerColor(kBlack);
	}
	else if(planenum == 2)
	{
		WaveformHist1->SetLineColor(kRed);
		WaveformHist1->SetMarkerColor(kRed);
		WaveformHist2->SetLineColor(kBlack);
		WaveformHist2->SetMarkerColor(kBlack);
	}
	
	gStyle->SetOptStat(0);
	gStyle->SetTitleSize(0.055,"T");

	TCanvas c1;
	c1.cd();
	WaveformHist1->Draw("HIST");
	WaveformHist2->Draw("HISTsame");
	WaveformHist1->SetTitle(Form("Plane %d, MC East vs West TPC: Average Waveform @ Anode",planenum));
	//WaveformHist1->SetTitle(Form("Plane %d, %s TPC: Average Waveform @ Anode",planenum,runtype1));
	WaveformHist1->GetXaxis()->SetTitle("Time Offset [#mus]");
	WaveformHist1->GetXaxis()->SetTitleSize(0.045);
	WaveformHist1->GetXaxis()->SetTitleOffset(1.05);
	WaveformHist1->GetXaxis()->SetLabelSize(0.04);
	WaveformHist1->GetXaxis()->SetRangeUser(timeTickSF*(-1*max_ticks_plot-0.5),timeTickSF*(max_ticks_plot+0.5));
	if(planenum == 0)
	{
		WaveformHist1->GetYaxis()->SetRangeUser(1.15*min(WaveformHist1->GetMinimum(),WaveformHist2->GetMinimum()),max(-0.5*min(WaveformHist1->GetMinimum(),WaveformHist2->GetMinimum()),1.2*max(WaveformHist1->GetMaximum(),WaveformHist2->GetMaximum())));
	}
	else if(planenum == 1)
	{
		WaveformHist1->GetYaxis()->SetRangeUser(1.2*min(WaveformHist1->GetMinimum(),WaveformHist2->GetMinimum()),1.2*max(WaveformHist1->GetMaximum(),WaveformHist2->GetMaximum()));
	}
	else if(planenum == 2)
	{
		WaveformHist1->GetYaxis()->SetRangeUser(-0.1*max(WaveformHist1->GetMaximum(),WaveformHist2->GetMaximum()),1.1*max(WaveformHist1->GetMaximum(),WaveformHist2->GetMaximum()));
	}
	WaveformHist1->GetYaxis()->SetTitle("Arb. Units");
	WaveformHist1->GetYaxis()->SetTitleSize(0.045);
	WaveformHist1->GetYaxis()->SetTitleOffset(1.12);
	WaveformHist1->GetYaxis()->SetLabelSize(0.04);
	TLegend *leg = new TLegend(0.60,0.72,0.85,0.85);
	leg->SetLineColor(kWhite);
	leg->SetTextSize(0.043);
	leg->AddEntry(WaveformHist1, "East TPC","L");
	leg->AddEntry(WaveformHist2, "West TPC","L");
	leg->Draw("SAME");

	//c1.SaveAs(Form("CombPlotWaveform_Anode_%s_Plane%d_%s_%dto%d.png",savefiletext,planenum,runtype1,angle_low,angle_high));
	c1.SaveAs(Form("Anode_%s_Plane%d_%dto%d_data.png",savefiletext,planenum,angle_low,angle_high));
}
