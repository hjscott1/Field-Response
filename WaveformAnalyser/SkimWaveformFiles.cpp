#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>

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

//DEFINES THE OUTPUT DIRECTORY
//const char* outputdir = (const char*) "/icarus/data/users/mrmooney/WaveformStudyResults/skimmed";
const char* outputdir = (const char*) "./";

//DETECTOR PARAMETERS
const int timePadding = 200;
const vector<float> max_nonlin_anode = {1.2, 0.9, 1.2};
const vector<float> max_nonlin_cathode = {1.2, 0.9, 1.2};

//DEFINES AddFiles METHODS
void AddFiles(TChain *ch, const char *dir, const char* substr);
void AddFilesList(TChain *ch, const char *listname, const char* substr);

//MAIN
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
	char *filetext;

	//DEFINE INPUT PARAMETERS
	if(argc < 2)
	{
		cout << endl << "Not enough input parameters provided.  Aborting." << endl << endl;
		return -1;
	}
	else
	{
		filepath = (char*) argv[1];
		if(argc == 3)
		{
			filetext = (char*) argv[2];
		}
		else
		{
			filetext = (char*) "root";
		}
	}

	bool useinputfilelist = false;
	TString teststring = (TString) filepath;
	if(teststring.Contains(".txt"))
	{
		useinputfilelist = true;
	}

	TChain* inputfiles = new TChain("trackTree");
	if(useinputfilelist == false)
	{
		AddFiles(inputfiles,filepath,filetext);
	}
	else
	{
		AddFilesList(inputfiles,filepath,filetext);
	}

	//CREATE OUTPUT FILE
	TFile outfile(Form("%s/wfstudy_skim.root",outputdir),"RECREATE");
	outfile.cd();

	int run, event, cryo, tpc, plane;
	float phi, theta;
	float anode_T, anode_W, anode_X, anode_Y, anode_Z, anode_NL;
	vector<short> anode_WFvals;
	float cathode_T, cathode_W, cathode_X, cathode_Y, cathode_Z, cathode_NL;
	vector<short> cathode_WFvals;

	float drift_T;

	//LOOK THROUGH ALL BRANCHES OF wfstudy.root
	TTree *trackTree = new TTree("trackTree","");
	trackTree->Branch("run",&run);
	trackTree->Branch("event",&event);
	trackTree->Branch("cryo",&cryo);
	trackTree->Branch("tpc",&tpc);
	trackTree->Branch("plane",&plane);
	trackTree->Branch("phi",&phi);
	trackTree->Branch("theta",&theta);
	trackTree->Branch("anode_T",&anode_T);
	trackTree->Branch("anode_W",&anode_W);
	trackTree->Branch("anode_X",&anode_X);
	trackTree->Branch("anode_Y",&anode_Y);
	trackTree->Branch("anode_Z",&anode_Z);
	trackTree->Branch("anode_NL",&anode_NL);
	trackTree->Branch("anode_WFvals",&anode_WFvals);
	trackTree->Branch("cathode_T",&cathode_T);
	trackTree->Branch("cathode_W",&cathode_W);
	trackTree->Branch("cathode_X",&cathode_X);
	trackTree->Branch("cathode_Y",&cathode_Y);
	trackTree->Branch("cathode_Z",&cathode_Z);
	trackTree->Branch("cathode_NL",&cathode_NL);
	trackTree->Branch("cathode_WFvals",&cathode_WFvals);
	trackTree->Branch("drift_T",&drift_T);

	for(int k = 0; k < 11*(2*timePadding+1); k++)
	{
		anode_WFvals.push_back(0);
		cathode_WFvals.push_back(0);
	}

	TTreeReader readerTracks(inputfiles);
	TTreeReaderValue<int> input_run(readerTracks, "run");
	TTreeReaderValue<int> input_event(readerTracks, "event");
	TTreeReaderValue<int> input_cryo(readerTracks, "cryo");
	TTreeReaderValue<int> input_tpc(readerTracks, "tpc");
	TTreeReaderValue<int> input_plane(readerTracks, "plane");
	TTreeReaderValue<float> input_phi(readerTracks, "phi");
	TTreeReaderValue<float> input_theta(readerTracks, "theta");
	TTreeReaderValue<float> input_anode_T(readerTracks, "anode_T");
	TTreeReaderValue<float> input_anode_W(readerTracks, "anode_W");
	TTreeReaderValue<float> input_anode_X(readerTracks, "anode_X");
	TTreeReaderValue<float> input_anode_Y(readerTracks, "anode_Y");
	TTreeReaderValue<float> input_anode_Z(readerTracks, "anode_Z");
	TTreeReaderValue<float> input_anode_NL(readerTracks, "anode_NL");
	TTreeReaderArray<short> input_anode_WFvals(readerTracks, "anode_WFvals");
	TTreeReaderValue<float> input_cathode_T(readerTracks, "cathode_T");
	TTreeReaderValue<float> input_cathode_W(readerTracks, "cathode_W");
	TTreeReaderValue<float> input_cathode_X(readerTracks, "cathode_X");
	TTreeReaderValue<float> input_cathode_Y(readerTracks, "cathode_Y");
	TTreeReaderValue<float> input_cathode_Z(readerTracks, "cathode_Z");
	TTreeReaderValue<float> input_cathode_NL(readerTracks, "cathode_NL");
	TTreeReaderArray<short> input_cathode_WFvals(readerTracks, "cathode_WFvals");

	int selected_tracks = 0;
	unsigned long long current_entry = 1;
	unsigned long long total_entries = readerTracks.GetEntries(true);
	cout << "Total Number of Input Tracks:  " << total_entries << endl;

	//DO SOME CUTS I GUESS...?
	while(readerTracks.Next())
	{
		if(current_entry % 1000 == 0)
		{
			cout << current_entry << "/" << total_entries << endl;
		}
		current_entry++;

		drift_T = *input_cathode_T - *input_anode_T;

		if(((*input_plane==0 && *input_anode_NL>max_nonlin_anode[0]) || (*input_plane==1 && *input_anode_NL>max_nonlin_anode[1]) || (*input_plane==2 && *input_anode_NL>max_nonlin_anode[2])) && ((*input_plane==0 && *input_cathode_NL>max_nonlin_cathode[0]) || (*input_plane==1 && *input_cathode_NL>max_nonlin_cathode[1]) || (*input_plane==2 && *input_cathode_NL>max_nonlin_cathode[2]))) continue;

		for(int k = 0; k < 11*(2*timePadding+1); k++)
		{
			anode_WFvals[k] = input_anode_WFvals[k];
			cathode_WFvals[k] = input_cathode_WFvals[k];
		}

		run = *input_run;
		event = *input_event;
		cryo = *input_cryo;
		tpc = *input_tpc;
		plane = *input_plane;

		phi = *input_phi;
		theta = *input_theta;

		anode_T = *input_anode_T;
		anode_W = *input_anode_W;
		anode_X = *input_anode_X;
		anode_Y = *input_anode_Y;
		anode_Z = *input_anode_Z;
		anode_NL = *input_anode_NL;

		cathode_T = *input_cathode_T;
		cathode_W = *input_cathode_W;
		cathode_X = *input_cathode_X;
		cathode_Y = *input_cathode_Y;
		cathode_Z = *input_cathode_Z;
		cathode_NL = *input_cathode_NL;

		trackTree->Fill();
		selected_tracks++;
	}

	cout << "Total Number of Skimmed Tracks:  " << selected_tracks << endl;

	outfile.cd();
	trackTree->Write();
	outfile.Close();

	timer.Stop();
	cout << "SkimWaveformFiles Runtime:  " << timer.RealTime() << " sec." << endl;

	return 0;
}

void AddFiles(TChain *ch, const char *dir, const char* substr)
{
	TString dirname = (TString) dir;
	dirname.ReplaceAll("/pnfs/","root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/");

	TSystemDirectory thisdir(dir, dir);
	TList *files = thisdir.GetListOfFiles();

	if(files)
	{
		TSystemFile *file;
		TString fname;
		TIter next(files);

		int counter = 1;
		while((file = (TSystemFile*)next()))
		{
			fname = file->GetName();
			if(!file->IsDirectory() && fname.Contains(substr) && fname.EndsWith(".root"))
			{
				ch->AddFile(Form("%s/%s",dirname.Data(),fname.Data()));
				counter++;
			}
		}
	}

	return;
}

void AddFilesList(TChain *ch, const char *listname, const char* substr)
{
	ifstream infile(listname);
	string filename;

	int counter = 1;
	while(getline(infile, filename))
	{
		TString fname = (TString) filename;
		if(fname.Contains(substr) && fname.EndsWith(".root"))
		{
			ch->AddFile(filename.c_str());
			counter++;
		}
	}

	return;
}
