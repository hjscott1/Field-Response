#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

#include "TROOT.h"
#include "ROOT/TThreadedObject.hxx"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TCut.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
#include "TList.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVirtualFFT.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystemDirectory.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

const char* mapfilepath = "map_yz_dqdx.root";
const int timePadding = 200;
const vector<float> max_nonlin_anode = {4.0, 4.0, 4.0}; // 1.2, 0.9, 1.2
const vector<float> max_nonlin_cathode = {4.0, 4.0, 4.0}; // 1.2, 0.9, 1.2

const float angle_min = 20.0;
const float angle_max = 80.0;
const float angle_step = 2.0;

const bool doBaselineCorr = true;
const bool trimExtrema = true;

const float timeTickSF = 0.5;

vector<short> TrimVecExtrema(vector<short> inputvec);

int main(int argc, char** argv)
{
	TStopwatch timer;
	timer.Start();

	gStyle->SetTitleSize(0.065,"T");
	//gErrorIgnoreLevel = kError;
	gErrorIgnoreLevel = kFatal;
	//gErrorIgnoreLevel = kWarning;
	double stops[5] = {0.00,0.34,0.61,0.84,1.00};
	double red[5] = {0.00,0.00,0.87,1.00,0.51};
	double green[5] = {0.00,0.81,1.00,0.20,0.00};
	double blue[5] = {0.51,1.00,0.12,0.00,0.00};
	TColor::CreateGradientColorTable(5,stops,red,green,blue,255);
	gStyle->SetNumberContours(255);
	gStyle->SetOptStat(0);
	TH1::AddDirectory(kFALSE);

	char *filepath;
	char *runtype;
	int planenum;
	int mapbins;
	int baselineMode; // 0 = use both sides, 1 = use left side only, 2 = use right side only

	if(argc < 6)
	{
		cout << endl << "Not enough input parameters provided.  Aborting." << endl << endl;
		return -1;
	}
	else
	{
		filepath = (char*) argv[1];
		runtype = (char*) argv[2];
		planenum = atoi(argv[3]);
		mapbins = atoi(argv[4]);
		baselineMode = atoi(argv[5]);
	}

	if((mapbins != 1) && (mapbins != 5) && (mapbins != 10) && (mapbins != 15))
	{
		cout << endl << "Unsupported number of bins for input map file.  Aborting." << endl << endl;
		return -2;
	}

	vector<pair<float, float>> phi_bins((int) (angle_max-angle_min)/angle_step);
	generate(phi_bins.begin(), phi_bins.end(), [] { static float x = angle_min-angle_step; float y = x+2*angle_step; return make_pair(x+=angle_step, y);} );


	std::string subfilepath(filepath);
	std::string runformat;
	
	std::size_t found_mc = subfilepath.find("wfstudy_mc");
	std::size_t found_data = subfilepath.find("wfstudy_data");
	if (found_mc!=std::string::npos)
	{
		runformat = "mc";
        cout << "Run format = " << runformat << endl;
	}
	else if (found_data!=std::string::npos)
	{
		runformat = "data";
        cout << "Run format = " << runformat << endl;
	}
	else
	{
		cout << endl << "Input file does not fit \"wfstudy_mc/data.root\" format. Aborting." << endl << endl;
		return -3;
	}

	size_t pos = subfilepath.find_last_of("/");
    std::string filedir = subfilepath.erase(pos);
	cout << "filedir = " << filedir << endl;

	std::string outdir = Form("%s/WFresults_Plane%i_%s_%s.root",filedir.c_str(),planenum,runtype,runformat.c_str());
	cout << "outdir = " << outdir << endl;

	TFile outfile(outdir.c_str(), "RECREATE");
	outfile.cd();

	TH1F* PhiHist1D = new TH1F("PhiHist1D","",91,-0.5,90.5);
	TH1F* ThetaHist1D = new TH1F("ThetaHist1D","",91,-0.5,90.5);

	TH1F* AnodeNonlinearityHist1D = new TH1F("AnodeNonlinearityHist1D","",100,0.0,5.0);
	TH2F* AnodeNonlinearityPhiHist2D = new TH2F("AnodeNonlinearityPhiHist2D","",91,-0.5,90.5,100,0.0,5.0);
	TH1F* CathodeNonlinearityHist1D = new TH1F("CathodeNonlinearityHist1D","",100,0.0,5.0);
	TH2F* CathodeNonlinearityPhiHist2D = new TH2F("CathodeNonlinearityPhiHist2D","",91,-0.5,90.5,100,0.0,5.0);

	TH1F* AnodeRecoHist1D[phi_bins.size()][mapbins];
	TH1F* CathodeRecoHist1D[phi_bins.size()][mapbins];
	TH2F* AnodeRecoHist2D[phi_bins.size()][mapbins];
	TH2F* CathodeRecoHist2D[phi_bins.size()][mapbins];
	TH3F* AnodeRecoHist3D[phi_bins.size()][mapbins];
	TH3F* CathodeRecoHist3D[phi_bins.size()][mapbins];

	for(int i = 0; i < phi_bins.size(); i++)
	{
		for(int j = 0; j < mapbins; j++)
		{
			if(mapbins != 1)
			{
				AnodeRecoHist1D[i][j] = new TH1F(Form("AnodeRecoHist1D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				CathodeRecoHist1D[i][j] = new TH1F(Form("CathodeRecoHist1D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));  
				AnodeRecoHist2D[i][j] = new TH2F(Form("AnodeRecoHist2D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
				CathodeRecoHist2D[i][j] = new TH2F(Form("CathodeRecoHist2D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
				AnodeRecoHist3D[i][j] = new TH3F(Form("AnodeRecoHist3D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
				CathodeRecoHist3D[i][j] = new TH3F(Form("CathodeRecoHist3D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
			}
			else
			{
				AnodeRecoHist1D[i][j] = new TH1F(Form("AnodeRecoHist1D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				CathodeRecoHist1D[i][j] = new TH1F(Form("CathodeRecoHist1D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));  
				AnodeRecoHist2D[i][j] = new TH2F(Form("AnodeRecoHist2D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
				CathodeRecoHist2D[i][j] = new TH2F(Form("CathodeRecoHist2D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
				AnodeRecoHist3D[i][j] = new TH3F(Form("AnodeRecoHist3D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
				CathodeRecoHist3D[i][j] = new TH3F(Form("CathodeRecoHist3D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),100,0.5,100.5);
			}
		}
	}

	TH2F* AnodeTrackHist2D[phi_bins.size()][mapbins];
	TH2F* CathodeTrackHist2D[phi_bins.size()][mapbins];
	TH2F* AnodeTrackUncertHist2D[phi_bins.size()][mapbins];
	TH2F* CathodeTrackUncertHist2D[phi_bins.size()][mapbins];
	TH3F* AnodeTrackHist3D[phi_bins.size()][mapbins];
	TH3F* CathodeTrackHist3D[phi_bins.size()][mapbins];

	for(int i = 0; i < phi_bins.size(); i++)
	{
		for(int j = 0; j < mapbins; j++)
		{
			if(mapbins != 1)
			{
				AnodeTrackHist2D[i][j] = new TH2F(Form("AnodeTrackHist2D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				CathodeTrackHist2D[i][j] = new TH2F(Form("CathodeTrackHist2D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				AnodeTrackUncertHist2D[i][j] = new TH2F(Form("AnodeTrackUncertHist2D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				CathodeTrackUncertHist2D[i][j] = new TH2F(Form("CathodeTrackUncertHist2D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				AnodeTrackHist3D[i][j] = new TH3F(Form("AnodeTrackHist3D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),101,-51,151);
				CathodeTrackHist3D[i][j] = new TH3F(Form("CathodeTrackHist3D_%.0fto%.0f_bin%d",phi_bins[i].first,phi_bins[i].second,j+1),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),101,-51,151);
			}
			else
			{
				AnodeTrackHist2D[i][j] = new TH2F(Form("AnodeTrackHist2D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				CathodeTrackHist2D[i][j] = new TH2F(Form("CathodeTrackHist2D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				AnodeTrackUncertHist2D[i][j] = new TH2F(Form("AnodeTrackUncertHist2D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				CathodeTrackUncertHist2D[i][j] = new TH2F(Form("CathodeTrackUncertHist2D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5));
				AnodeTrackHist3D[i][j] = new TH3F(Form("AnodeTrackHist3D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),101,-51,151);
				CathodeTrackHist3D[i][j] = new TH3F(Form("CathodeTrackHist3D_%.0fto%.0f",phi_bins[i].first,phi_bins[i].second),"",11,-5.5,5.5,2*timePadding+1,timeTickSF*(-1*timePadding-0.5),timeTickSF*(timePadding+0.5),101,-51,151);
			}
		}
	}

	int tpcNum;
	if(!strcmp(runtype,"X"))
	{
		tpcNum = -1;
	}
	else if(!strcmp(runtype,"E"))
	{
		tpcNum = 0;
	}
	else if(!strcmp(runtype,"W"))
	{
		tpcNum = 1;
	}

	TChain* inputfile = new TChain("trackTree");
	inputfile->Add(filepath);

	TTreeReader readerTracks(inputfile);
	TTreeReaderValue<int> run(readerTracks, "run");
	TTreeReaderValue<int> event(readerTracks, "event");
	TTreeReaderValue<int> cryo(readerTracks, "cryo");
	TTreeReaderValue<int> tpc(readerTracks, "tpc");
	TTreeReaderValue<int> plane(readerTracks, "plane");
	TTreeReaderValue<float> phi(readerTracks, "phi");
	TTreeReaderValue<float> theta(readerTracks, "theta");
	TTreeReaderValue<float> anode_T(readerTracks, "anode_T");
	TTreeReaderValue<float> anode_W(readerTracks, "anode_W");
	TTreeReaderValue<float> anode_X(readerTracks, "anode_X");
	TTreeReaderValue<float> anode_Y(readerTracks, "anode_Y");
	TTreeReaderValue<float> anode_Z(readerTracks, "anode_Z");
	TTreeReaderValue<float> anode_NL(readerTracks, "anode_NL");
	TTreeReaderArray<short> anode_WFvals(readerTracks, "anode_WFvals");
	TTreeReaderValue<float> cathode_T(readerTracks, "cathode_T");
	TTreeReaderValue<float> cathode_W(readerTracks, "cathode_W");
	TTreeReaderValue<float> cathode_X(readerTracks, "cathode_X");
	TTreeReaderValue<float> cathode_Y(readerTracks, "cathode_Y");
	TTreeReaderValue<float> cathode_Z(readerTracks, "cathode_Z");
	TTreeReaderValue<float> cathode_NL(readerTracks, "cathode_NL");
	TTreeReaderArray<short> cathode_WFvals(readerTracks, "cathode_WFvals");

	TFile* mapfile = new TFile(mapfilepath,"READ");

	//TH2F* ScaleMapRaw[3][2];
	//TH2F* ScaleMapBins[3][2];
	TH2F* ScaleMapRaw[3][4];
    TH2F* ScaleMapBins[3][4];
	
	//vector<string> tpc_names = {"E", "W"};
	vector<string> tpc_names = {"EE", "EW", "WE", "WW"};

	for(int i = 0; i < 3; i++)
	{
		//for(int j = 0; j < 2; j++)
		for(int j = 0; j < 4; j++)
		{
			ScaleMapRaw[i][j] = (TH2F*) mapfile->Get(Form("map_yz_dqdx_plane%d_%s",i,tpc_names[j].c_str()));
			if(mapbins == 5)
			{
				ScaleMapBins[i][j] = (TH2F*) mapfile->Get(Form("map_yz_dqdx_5bin_plane%d_%s",i,tpc_names[j].c_str()));
			}
			else if(mapbins == 10)
			{
				ScaleMapBins[i][j] = (TH2F*) mapfile->Get(Form("map_yz_dqdx_10bin_plane%d_%s",i,tpc_names[j].c_str()));
			}
			else
			{
				ScaleMapBins[i][j] = (TH2F*) mapfile->Get(Form("map_yz_dqdx_15bin_plane%d_%s",i,tpc_names[j].c_str()));
			}

			if(mapbins == 1)
			{
				for(int h = 1; h <= ScaleMapBins[i][j]->GetNbinsX(); h++)
				{
					for(int k = 1; k <= ScaleMapBins[i][j]->GetNbinsY(); k++)
					{
						ScaleMapBins[i][j]->SetBinContent(h,k,1);
					}
				}
			}
		}
	}

	int counter_A[phi_bins.size()][mapbins] = {0};
	int counter_C[phi_bins.size()][mapbins] = {0};

	auto all_anode_WFvals = vector<vector<vector<vector<vector<short>>>>>(phi_bins.size(),vector<vector<vector<vector<short>>>>(mapbins,vector<vector<vector<short>>>(11,vector<vector<short>>(2*timePadding+1,vector<short>(0)))));
	auto all_cathode_WFvals = vector<vector<vector<vector<vector<short>>>>>(phi_bins.size(),vector<vector<vector<vector<short>>>>(mapbins,vector<vector<vector<short>>>(11,vector<vector<short>>(2*timePadding+1,vector<short>(0)))));

	unsigned long long current_entry = 1;
	unsigned long long total_entries = readerTracks.GetEntries(true);
	cout << "Total Number of Input Tracks:  " << total_entries << endl;

	while(readerTracks.Next())
	{
		if(current_entry % 10000 == 0)
		{
			cout << current_entry << "/" << total_entries << endl;
		}
		current_entry++;

		if(*plane != planenum || ((tpcNum != -1) && (*tpc != tpcNum))) continue;
		
		int phi_index = -1;
		for(int i = 0; i < phi_bins.size(); i++)
		{
			if((*phi > phi_bins[i].first) && (*phi < phi_bins[i].second))
			{
				phi_index = i;
			}
		}

		if(phi_index < 0)
		{
			continue;
		}

		int bin_anode = round(ScaleMapBins[*plane][2*(*cryo)+(*tpc)]->Interpolate(*anode_Z,*anode_Y));
		int bin_cathode = round(ScaleMapBins[*plane][2*(*cryo)+(*tpc)]->Interpolate(*cathode_Z,*cathode_Y));

		if(bin_anode < 1)
		{
			bin_anode = 1;
		}
		if(bin_cathode < 1)
		{
			bin_cathode = 1;
		}
		if(bin_anode > mapbins)
		{
			bin_anode = mapbins;
		}
		if(bin_cathode > mapbins)
		{
			bin_cathode = mapbins;
		}

		if(*anode_NL < max_nonlin_anode[planenum])
		{
			for(int q = 0; q < 11; q++)
			{
				for(int k = 0; k < 2*timePadding+1; k++)
				{
					all_anode_WFvals[phi_index][bin_anode-1][q][k].push_back(anode_WFvals[q*(2*timePadding+1)+k]);

					if((counter_A[phi_index][bin_anode-1] < 100) && (q == 5))
					{
						AnodeRecoHist2D[phi_index][bin_anode-1]->SetBinContent(k+1,counter_A[phi_index][bin_anode-1]+1,anode_WFvals[q*(2*timePadding+1)+k]);
					}
					if(counter_A[phi_index][bin_anode-1] < 100)
					{
						AnodeRecoHist3D[phi_index][bin_anode-1]->SetBinContent(q+1,k+1,counter_A[phi_index][bin_anode-1]+1,anode_WFvals[q*(2*timePadding+1)+k]);
					}
				}
			}

			PhiHist1D->Fill(*phi);
			ThetaHist1D->Fill(*theta);
			AnodeNonlinearityHist1D->Fill(*anode_NL);
			AnodeNonlinearityPhiHist2D->Fill(*phi,*anode_NL);

			counter_A[phi_index][bin_anode-1]++;
		}

		if(*cathode_NL < max_nonlin_cathode[planenum])
		{
			for(int q = 0; q < 11; q++)
			{
				for(int k = 0; k < 2*timePadding+1; k++)
				{
					all_cathode_WFvals[phi_index][bin_cathode-1][q][k].push_back(cathode_WFvals[q*(2*timePadding+1)+k]);

					if((counter_C[phi_index][bin_cathode-1] < 100) && (q == 5))
					{
						CathodeRecoHist2D[phi_index][bin_cathode-1]->SetBinContent(k+1,counter_C[phi_index][bin_cathode-1]+1,cathode_WFvals[q*(2*timePadding+1)+k]);
					}
					if(counter_C[phi_index][bin_cathode-1] < 100)
					{
						CathodeRecoHist3D[phi_index][bin_cathode-1]->SetBinContent(q+1,k+1,counter_C[phi_index][bin_cathode-1]+1,cathode_WFvals[q*(2*timePadding+1)+k]);
					}
				}
			}

			PhiHist1D->Fill(*phi);
			ThetaHist1D->Fill(*theta);
			CathodeNonlinearityHist1D->Fill(*cathode_NL);
			CathodeNonlinearityPhiHist2D->Fill(*phi,*cathode_NL);

			counter_C[phi_index][bin_cathode-1]++;
		}
	}

	int total_tracks_A = 0;
	int total_tracks_C = 0;
	for(int i = 0; i < phi_bins.size(); i++)
	{
		for(int j = 0; j < mapbins; j++)
		{
			total_tracks_A += all_anode_WFvals[i][j][5][0].size();
			total_tracks_C += all_cathode_WFvals[i][j][5][0].size();
		}
	}

	cout << "Plane " << planenum << " Anode Tracks:  " << total_tracks_A << endl;
	cout << "Plane " << planenum << " Cathode Tracks:  " << total_tracks_C << endl;

	for(int i = 0; i < phi_bins.size(); i++)
	{
		for(int j = 0; j < mapbins; j++)
		{
			for(int q = 0; q < 11; q++)
			{
				for(int k = 0; k < 2*timePadding+1; k++)
				{
					vector<short> temp_all_anode_WFvals;
					vector<short> temp_all_cathode_WFvals;
					if(trimExtrema == true)
					{
						temp_all_anode_WFvals = TrimVecExtrema(all_anode_WFvals[i][j][q][k]);
						temp_all_cathode_WFvals = TrimVecExtrema(all_cathode_WFvals[i][j][q][k]);
					}
					else
					{
						temp_all_anode_WFvals = all_anode_WFvals[i][j][q][k];
						temp_all_cathode_WFvals = all_cathode_WFvals[i][j][q][k];
					}

					float anodeSum = 0.0;
					float cathodeSum = 0.0;
					float anodeSumSq = 0.0;
					float cathodeSumSq = 0.0;
					float anodeN = 0.0;
					float cathodeN = 0.0;

					const int numEntriesA = temp_all_anode_WFvals.size();
					for(int h = 0; h < numEntriesA; h++)
					{
						AnodeTrackHist3D[i][j]->Fill(q-5,timeTickSF*(k-timePadding),temp_all_anode_WFvals[h]);
						anodeSum += temp_all_anode_WFvals[h];
						anodeSumSq += pow(temp_all_anode_WFvals[h],2);
						anodeN += 1.0;
					}

					const int numEntriesC = temp_all_cathode_WFvals.size();
					for(int h = 0; h < numEntriesC; h++)
					{
						CathodeTrackHist3D[i][j]->Fill(q-5,timeTickSF*(k-timePadding),temp_all_cathode_WFvals[h]);
						cathodeSum += temp_all_cathode_WFvals[h];
						cathodeSumSq += pow(temp_all_cathode_WFvals[h],2);
						cathodeN += 1.0;
					}

					if(anodeN != 0)
					{
						AnodeTrackHist2D[i][j]->SetBinContent(q+1,k+1,anodeSum/anodeN);
						AnodeTrackUncertHist2D[i][j]->SetBinContent(q+1,k+1,sqrt(anodeSumSq/anodeN - pow(anodeSum/anodeN,2))/sqrt(anodeN));
					}
					else
					{
						AnodeTrackHist2D[i][j]->SetBinContent(q+1,k+1,0.0);
						AnodeTrackUncertHist2D[i][j]->SetBinContent(q+1,k+1,0.0);
					}

					if(cathodeN != 0)
					{
						CathodeTrackHist2D[i][j]->SetBinContent(q+1,k+1,cathodeSum/cathodeN);
						CathodeTrackUncertHist2D[i][j]->SetBinContent(q+1,k+1,sqrt(cathodeSumSq/cathodeN - pow(cathodeSum/cathodeN,2))/sqrt(cathodeN));
					}
					else
					{
						CathodeTrackHist2D[i][j]->SetBinContent(q+1,k+1,0.0);
						CathodeTrackUncertHist2D[i][j]->SetBinContent(q+1,k+1,0.0);
					}
				}
			}

			for(int k = 0; k < 2*timePadding+1; k++)
			{
				vector<short> temp_all_anode_WFvals;
				vector<short> temp_all_cathode_WFvals;
				if(trimExtrema == true)
				{
					temp_all_anode_WFvals = TrimVecExtrema(all_anode_WFvals[i][j][5][k]);
					temp_all_cathode_WFvals = TrimVecExtrema(all_cathode_WFvals[i][j][5][k]);
				}
				else
				{
					temp_all_anode_WFvals = all_anode_WFvals[i][j][5][k];
					temp_all_cathode_WFvals = all_cathode_WFvals[i][j][5][k];
				}

				float anodeSum = 0.0;
				float cathodeSum = 0.0; 
				float anodeN = 0.0;
				float cathodeN = 0.0; 

				const int numEntriesA = temp_all_anode_WFvals.size();
				for(int h = 0; h < numEntriesA; h++)
				{
					anodeSum += temp_all_anode_WFvals[h];
					anodeN += 1.0;
				}

				const int numEntriesC = temp_all_cathode_WFvals.size();
				for(int h = 0; h < numEntriesC; h++)
				{
					cathodeSum += temp_all_cathode_WFvals[h];
					cathodeN += 1.0;
				}

				if(anodeN != 0)
				{
					AnodeRecoHist1D[i][j]->SetBinContent(k+1,anodeSum/anodeN);
				}
				else
				{
					AnodeRecoHist1D[i][j]->SetBinContent(k+1,0.0);
				}

				if(cathodeN != 0)
				{
					CathodeRecoHist1D[i][j]->SetBinContent(k+1,cathodeSum/cathodeN);
				}
				else
				{
					CathodeRecoHist1D[i][j]->SetBinContent(k+1,0.0);
				}
			}
		}
	}

	if(doBaselineCorr == true)
	{
		for(int i = 0; i < phi_bins.size(); i++)
		{
			for(int j = 0; j < mapbins; j++)
			{
				TH1F *AnodeTrackHist2D_projY = (TH1F*) AnodeTrackHist2D[i][j]->ProjectionY();
				TH1F *CathodeTrackHist2D_projY = (TH1F*) CathodeTrackHist2D[i][j]->ProjectionY();

				float anodeCorr_left_avgTime = 0.0;
				float anodeCorr_left_avgVal = 0.0;
				float anodeCorr_left_N = 0.0;
				float anodeCorr_right_avgTime = 0.0;
				float anodeCorr_right_avgVal = 0.0;
				float anodeCorr_right_N = 0.0;

				float cathodeCorr_left_avgTime = 0.0;
				float cathodeCorr_left_avgVal = 0.0;
				float cathodeCorr_left_N = 0.0;
				float cathodeCorr_right_avgTime = 0.0;
				float cathodeCorr_right_avgVal = 0.0;
				float cathodeCorr_right_N = 0.0;

				for(int q = 0; q < 11; q++)
				{
					for(int k = 0; k < 2*timePadding+1; k++)
					{
						if((k < 1*timePadding/5) && (q > 8))
						{
							anodeCorr_left_avgTime += AnodeTrackHist2D_projY->GetBinCenter(k+1);
							anodeCorr_left_avgVal += AnodeTrackHist2D[i][j]->GetBinContent(q+1,k+1);
							anodeCorr_left_N += 1.0;

							cathodeCorr_left_avgTime += CathodeTrackHist2D_projY->GetBinCenter(k+1);
							cathodeCorr_left_avgVal += CathodeTrackHist2D[i][j]->GetBinContent(q+1,k+1);
							cathodeCorr_left_N += 1.0;
						}
						else if((k > 9*timePadding/5) && (q < 2))
						{
							anodeCorr_right_avgTime += AnodeTrackHist2D_projY->GetBinCenter(k+1);
							anodeCorr_right_avgVal += AnodeTrackHist2D[i][j]->GetBinContent(q+1,k+1);
							anodeCorr_right_N += 1.0;

							cathodeCorr_right_avgTime += CathodeTrackHist2D_projY->GetBinCenter(k+1);
							cathodeCorr_right_avgVal += CathodeTrackHist2D[i][j]->GetBinContent(q+1,k+1);
							cathodeCorr_right_N += 1.0;
						}
					}
				}

				anodeCorr_left_avgTime /= anodeCorr_left_N;
				anodeCorr_left_avgVal /= anodeCorr_left_N;
				anodeCorr_right_avgTime /= anodeCorr_right_N;
				anodeCorr_right_avgVal /= anodeCorr_right_N;

				cathodeCorr_left_avgTime /= cathodeCorr_left_N;
				cathodeCorr_left_avgVal /= cathodeCorr_left_N;
				cathodeCorr_right_avgTime /= cathodeCorr_right_N;
				cathodeCorr_right_avgVal /= cathodeCorr_right_N;

				if(baselineMode == 1)
				{
					anodeCorr_right_avgVal = anodeCorr_left_avgVal;
					cathodeCorr_right_avgVal = cathodeCorr_left_avgVal;
				}
				else if(baselineMode == 2)
				{
					anodeCorr_left_avgVal = anodeCorr_right_avgVal;
					cathodeCorr_left_avgVal = cathodeCorr_right_avgVal;
				}

				float baselineCorrAnode_slope = (anodeCorr_right_avgVal-anodeCorr_left_avgVal)/(anodeCorr_right_avgTime-anodeCorr_left_avgTime);
				float baselineCorrAnode_intercept = (anodeCorr_left_avgVal+anodeCorr_right_avgVal)/2.0;

				float baselineCorrCathode_slope = (cathodeCorr_right_avgVal-cathodeCorr_left_avgVal)/(cathodeCorr_right_avgTime-cathodeCorr_left_avgTime);
				float baselineCorrCathode_intercept = (cathodeCorr_left_avgVal+cathodeCorr_right_avgVal)/2.0;

				for(int q = 0; q < 11; q++)
				{
					for(int k = 0; k < 2*timePadding+1; k++)
					{
						AnodeTrackHist2D[i][j]->SetBinContent(q+1,k+1,AnodeTrackHist2D[i][j]->GetBinContent(q+1,k+1)-(baselineCorrAnode_slope*AnodeTrackHist2D_projY->GetBinCenter(k+1)+baselineCorrAnode_intercept));
						CathodeTrackHist2D[i][j]->SetBinContent(q+1,k+1,CathodeTrackHist2D[i][j]->GetBinContent(q+1,k+1)-(baselineCorrCathode_slope*CathodeTrackHist2D_projY->GetBinCenter(k+1)+baselineCorrCathode_intercept));
					}
				}

				anodeCorr_left_avgTime = 0.0;
				anodeCorr_left_avgVal = 0.0;
				anodeCorr_left_N = 0.0;
				anodeCorr_right_avgTime = 0.0;
				anodeCorr_right_avgVal = 0.0;
				anodeCorr_right_N = 0.0;

				cathodeCorr_left_avgTime = 0.0;
				cathodeCorr_left_avgVal = 0.0;
				cathodeCorr_left_N = 0.0;
				cathodeCorr_right_avgTime = 0.0;
				cathodeCorr_right_avgVal = 0.0;
				cathodeCorr_right_N = 0.0;

				for(int k = 0; k < 2*timePadding+1; k++)
				{
					if(k < 1*timePadding/5)
					{
						anodeCorr_left_avgTime += AnodeRecoHist1D[i][j]->GetBinCenter(k+1);
						cathodeCorr_left_avgTime += CathodeRecoHist1D[i][j]->GetBinCenter(k+1);
						anodeCorr_left_avgVal += AnodeRecoHist1D[i][j]->GetBinContent(k+1);
						cathodeCorr_left_avgVal += CathodeRecoHist1D[i][j]->GetBinContent(k+1);
						anodeCorr_left_N += 1.0;
						cathodeCorr_left_N += 1.0;
					}
					else if(k > 9*timePadding/5)
					{
						anodeCorr_right_avgTime += AnodeRecoHist1D[i][j]->GetBinCenter(k+1);
						cathodeCorr_right_avgTime += CathodeRecoHist1D[i][j]->GetBinCenter(k+1);
						anodeCorr_right_avgVal += AnodeRecoHist1D[i][j]->GetBinContent(k+1);
						cathodeCorr_right_avgVal += CathodeRecoHist1D[i][j]->GetBinContent(k+1);
						anodeCorr_right_N += 1.0;
						cathodeCorr_right_N += 1.0;
					}
				}

				anodeCorr_left_avgTime /= anodeCorr_left_N;
				cathodeCorr_left_avgTime /= cathodeCorr_left_N;
				anodeCorr_left_avgVal /= anodeCorr_left_N;
				cathodeCorr_left_avgVal /= cathodeCorr_left_N;
				anodeCorr_right_avgTime /= anodeCorr_right_N;
				cathodeCorr_right_avgTime /= cathodeCorr_right_N;
				anodeCorr_right_avgVal /= anodeCorr_right_N;
				cathodeCorr_right_avgVal /= cathodeCorr_right_N;

				if(baselineMode == 1)
				{
					anodeCorr_right_avgVal = anodeCorr_left_avgVal;
					cathodeCorr_right_avgVal = cathodeCorr_left_avgVal;
				}
				else if(baselineMode == 2)
				{
					anodeCorr_left_avgVal = anodeCorr_right_avgVal;
					cathodeCorr_left_avgVal = cathodeCorr_right_avgVal;
				}

				baselineCorrAnode_slope = (anodeCorr_right_avgVal-anodeCorr_left_avgVal)/(anodeCorr_right_avgTime-anodeCorr_left_avgTime);
				baselineCorrAnode_intercept = (anodeCorr_left_avgVal+anodeCorr_right_avgVal)/2.0;
				baselineCorrCathode_slope = (cathodeCorr_right_avgVal-cathodeCorr_left_avgVal)/(cathodeCorr_right_avgTime-cathodeCorr_left_avgTime);
				baselineCorrCathode_intercept = (cathodeCorr_left_avgVal+cathodeCorr_right_avgVal)/2.0;

				for(int k = 0; k < 2*timePadding+1; k++)
				{
					AnodeRecoHist1D[i][j]->SetBinContent(k+1,AnodeRecoHist1D[i][j]->GetBinContent(k+1)-(baselineCorrAnode_slope*AnodeRecoHist1D[i][j]->GetBinCenter(k+1)+baselineCorrAnode_intercept));
					CathodeRecoHist1D[i][j]->SetBinContent(k+1,CathodeRecoHist1D[i][j]->GetBinContent(k+1)-(baselineCorrCathode_slope*CathodeRecoHist1D[i][j]->GetBinCenter(k+1)+baselineCorrCathode_intercept));
				}
			}
		}
	}

	outfile.cd();

	PhiHist1D->Write();
	ThetaHist1D->Write();

	AnodeNonlinearityHist1D->Write();
	AnodeNonlinearityPhiHist2D->Write();
	CathodeNonlinearityHist1D->Write();
	CathodeNonlinearityPhiHist2D->Write();

	for(int i = 0; i < phi_bins.size(); i++)
	{
		for(int j = 0; j < mapbins; j++)
		{
			AnodeRecoHist1D[i][j]->Write();
			CathodeRecoHist1D[i][j]->Write();
			AnodeRecoHist2D[i][j]->Write();
			CathodeRecoHist2D[i][j]->Write();
			AnodeRecoHist3D[i][j]->Write();
			CathodeRecoHist3D[i][j]->Write();

			AnodeTrackHist2D[i][j]->Write();
			CathodeTrackHist2D[i][j]->Write();
			AnodeTrackUncertHist2D[i][j]->Write();
			CathodeTrackUncertHist2D[i][j]->Write();
			AnodeTrackHist3D[i][j]->Write();
			CathodeTrackHist3D[i][j]->Write();
		}
	}

	timer.Stop();
	cout << "ProcessWaveformFiles Runtime:  " << timer.RealTime() << " sec." << endl;

	return 0;
}

vector<short> TrimVecExtrema(vector<short> inputvec)
{
	sort(inputvec.begin(), inputvec.end());

	vector<short> newvec;
	for(int i = round(inputvec.size()/10.0); i < round(9.0*inputvec.size()/10.0); i++)
	{
		newvec.push_back(inputvec[i]);
	}

	return newvec;
}
