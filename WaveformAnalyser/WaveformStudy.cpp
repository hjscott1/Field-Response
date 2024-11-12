//ROOT INCLUDES
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <filesystem>

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

#include "TrackCaloSkimmerObj.h"

using namespace std;
//SET DETECTOR GEOMETRY AND PARAMETERS
const float minDistOffset = 2.0; //cm	old: 13.0
const float maxDistOffset = 5.0; //cm	old: 16.0
const int timePadding = 200; //ns
const vector<float> threshold = {6.0, 6.0, 6.0};
const int max_offset = 5;

const float driftVel = 0.1571; //Data
const float timeTickSF = 0.5; //mus
const float wirePitch = 0.3;
const int maxTicks = 3400; //ns

const float WF_top = 200; //cm
const float WF_bottom = -200; //cm
const float WF_upstream = 0; //cm
const float WF_downstream = 500; //cm
const float WF_cathode = 0; //cm
const float WF_ACdist = 200; //cm

//METHODS THAT ADD THE FILES TO THE SCRIPT USING INPUT PARAMETERS
void AddFiles(TChain *ch, const char *dir, const char* substr, int minfilenum, int maxfilenum);
void AddFilesList(TChain *ch, const char *listname, const char* substr, int minfilenum, int maxfilenum);
void CheckIfDirectoryExists(std::filesystem::path filepath);

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

	gSystem->Load("TrackCaloSkimmerObj_h.so");

	//INPUT PARAMETERS
	char *filepath;
	int filenums;
	char *filetext;
	char *runformat;
	char *outputfilename;

	//ASSIGNS THE INPUT PARAMETERS VARIABLES
	if(argc < 5)
	{
		cout << endl << "Not enough input parameters provided.  Aborting." << endl << endl;
		return -1;
	}
	else
	{
		filepath = (char*) argv[1];
		filenums = atoi(argv[2]);
		runformat = (char*) argv[3];
		outputfilename = (char*) argv[4];

		if(argc == 6)
		{
			filetext = (char*) argv[5];
		}
		else
		{
			filetext = (char*) "root";
		}
	}

	//ADDS THE FILES BASED ON THE INPUT PARAMETERS
	bool useinputfilelist = false;
	TString teststring = (TString) filepath;
	if(teststring.Contains(".list"))
	{
		useinputfilelist = true;
	}

	TChain* inputfiles = new TChain("caloskim/TrackCaloSkim");
	if(useinputfilelist == false)
	{
		AddFiles(inputfiles,filepath,filetext,1,filenums);
	}
	else
	{
		AddFilesList(inputfiles,filepath,filetext,1,filenums);
	}

	//OUTPUTS THE FILE
	std::string outputdir = Form("Results/%s/%s/wfstudy_%s.root",runformat,outputfilename,runformat);

	CheckIfDirectoryExists(outputdir);

	TFile outfile(outputdir.c_str(),"RECREATE");
	outfile.cd();

	//STORES THE TTREE BRANCHES AS PARAMETERS
	int run, event, cryo, tpc, plane;
	float phi, theta;
	float anode_T, anode_W, anode_X, anode_Y, anode_Z, anode_NL;
	vector<short> anode_WFvals;
	float cathode_T, cathode_W, cathode_X, cathode_Y, cathode_Z, cathode_NL;
	vector<short> cathode_WFvals;

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

	for(int k = 0; k < 11*(2*timePadding+1); k++)
	{
		anode_WFvals.push_back(0);
		cathode_WFvals.push_back(0);
	}

	TTreeReader readerTracks(inputfiles);
	TTreeReaderValue<int> runNum(readerTracks, "meta.run");
	TTreeReaderValue<int> eventNum(readerTracks, "meta.evt");
	TTreeReaderValue<int> tagType(readerTracks, "selected");
	TTreeReaderValue<int> cryoNum(readerTracks, "cryostat");
	TTreeReaderValue<float> trackDirX(readerTracks, "dir.x");
	TTreeReaderValue<float> trackDirY(readerTracks, "dir.y");
	TTreeReaderValue<float> trackDirZ(readerTracks, "dir.z");
	TTreeReaderValue<float> minTime_TPCE_plane0(readerTracks, "hit_min_time_p0_tpcE");
	TTreeReaderValue<float> maxTime_TPCE_plane0(readerTracks, "hit_max_time_p0_tpcE");
	TTreeReaderValue<float> minTime_TPCW_plane0(readerTracks, "hit_min_time_p0_tpcW");
	TTreeReaderValue<float> maxTime_TPCW_plane0(readerTracks, "hit_max_time_p0_tpcW");
	TTreeReaderArray<float> trackHitTimes_plane0(readerTracks, "hits0.h.time");
	TTreeReaderArray<unsigned short> trackHitWires_plane0(readerTracks, "hits0.h.wire");
	TTreeReaderArray<unsigned short> trackHitTPCs_plane0(readerTracks, "hits0.h.tpc");
	TTreeReaderArray<float> trackHitXvals_plane0(readerTracks, "hits0.tp.x");
	TTreeReaderArray<float> trackHitYvals_plane0(readerTracks, "hits0.tp.y");
	TTreeReaderArray<float> trackHitZvals_plane0(readerTracks, "hits0.tp.z");
	TTreeReaderArray<bool> trackHitIsOnTraj_plane0(readerTracks, "hits0.ontraj");
	TTreeReaderArray<sbn::WireInfo> trackWf_plane0(readerTracks, "wires0");
	TTreeReaderValue<float> minTime_TPCE_plane1(readerTracks, "hit_min_time_p1_tpcE");
	TTreeReaderValue<float> maxTime_TPCE_plane1(readerTracks, "hit_max_time_p1_tpcE");
	TTreeReaderValue<float> minTime_TPCW_plane1(readerTracks, "hit_min_time_p1_tpcW");
	TTreeReaderValue<float> maxTime_TPCW_plane1(readerTracks, "hit_max_time_p1_tpcW");
	TTreeReaderArray<float> trackHitTimes_plane1(readerTracks, "hits1.h.time");
	TTreeReaderArray<unsigned short> trackHitWires_plane1(readerTracks, "hits1.h.wire");
	TTreeReaderArray<unsigned short> trackHitTPCs_plane1(readerTracks, "hits1.h.tpc");
	TTreeReaderArray<float> trackHitXvals_plane1(readerTracks, "hits1.tp.x");
	TTreeReaderArray<float> trackHitYvals_plane1(readerTracks, "hits1.tp.y");
	TTreeReaderArray<float> trackHitZvals_plane1(readerTracks, "hits1.tp.z");
	TTreeReaderArray<bool> trackHitIsOnTraj_plane1(readerTracks, "hits1.ontraj");
	TTreeReaderArray<sbn::WireInfo> trackWf_plane1(readerTracks, "wires1");
	TTreeReaderValue<float> minTime_TPCE_plane2(readerTracks, "hit_min_time_p2_tpcE");
	TTreeReaderValue<float> maxTime_TPCE_plane2(readerTracks, "hit_max_time_p2_tpcE");
	TTreeReaderValue<float> minTime_TPCW_plane2(readerTracks, "hit_min_time_p2_tpcW");
	TTreeReaderValue<float> maxTime_TPCW_plane2(readerTracks, "hit_max_time_p2_tpcW");
	TTreeReaderArray<float> trackHitTimes_plane2(readerTracks, "hits2.h.time");
	TTreeReaderArray<unsigned short> trackHitWires_plane2(readerTracks, "hits2.h.wire");
	TTreeReaderArray<unsigned short> trackHitTPCs_plane2(readerTracks, "hits2.h.tpc");
	TTreeReaderArray<float> trackHitXvals_plane2(readerTracks, "hits2.tp.x");
	TTreeReaderArray<float> trackHitYvals_plane2(readerTracks, "hits2.tp.y");
	TTreeReaderArray<float> trackHitZvals_plane2(readerTracks, "hits2.tp.z");
	TTreeReaderArray<bool> trackHitIsOnTraj_plane2(readerTracks, "hits2.ontraj");
	TTreeReaderArray<sbn::WireInfo> trackWf_plane2(readerTracks, "wires2");

	float minTime_TPCE;
	float maxTime_TPCE;
	float minTime_TPCW;
	float maxTime_TPCW;
	vector<float> trackHitTimes;
	vector<unsigned short> trackHitWires;
	vector<unsigned short> trackHitTPCs;
	vector<float> trackHitXvals;
	vector<float> trackHitYvals;
	vector<float> trackHitZvals;
	vector<bool> trackHitIsOnTraj;
	vector<sbn::WireInfo> trackWf;

	int selected_tracks[3] = {0};

	unsigned long long current_entry = 1;
	unsigned long long total_entries = readerTracks.GetEntries(true);
	cout << "Total Number of Input Tracks:  " << total_entries << endl;
	
	while(readerTracks.Next())
	{	
		if(current_entry % 10000 == 0)
		{
			cout << "Event " << *eventNum << ". Hit: " << current_entry << "/" << total_entries << endl;
		}
		current_entry++;

		for(int planenum = 0; planenum <= 2; planenum++)
		{
			trackHitTimes.clear();
			trackHitWires.clear();
			trackHitTPCs.clear();
			trackHitXvals.clear();
			trackHitYvals.clear();
			trackHitZvals.clear();
			trackHitIsOnTraj.clear();
			trackWf.clear();

			if(planenum == 0)
			{
				minTime_TPCE = *minTime_TPCE_plane0;
				maxTime_TPCE = *maxTime_TPCE_plane0;
				minTime_TPCW = *minTime_TPCW_plane0;
				maxTime_TPCW = *maxTime_TPCW_plane0;
				for(int k = 0; k < trackHitTimes_plane0.GetSize(); k++)
				{
					trackHitTimes.push_back(trackHitTimes_plane0[k]);
					trackHitWires.push_back(trackHitWires_plane0[k]);
					trackHitTPCs.push_back(trackHitTPCs_plane0[k]);
					trackHitXvals.push_back(trackHitXvals_plane0[k]);
					trackHitYvals.push_back(trackHitYvals_plane0[k]);
					trackHitZvals.push_back(trackHitZvals_plane0[k]);
					trackHitIsOnTraj.push_back(trackHitIsOnTraj_plane0[k]);
				}
				for(int k = 0; k < trackWf_plane0.GetSize(); k++)
				{
					trackWf.push_back(trackWf_plane0[k]);
				}
			}
			else if(planenum == 1)
			{
				minTime_TPCE = *minTime_TPCE_plane1;
				maxTime_TPCE = *maxTime_TPCE_plane1;
				minTime_TPCW = *minTime_TPCW_plane1;
				maxTime_TPCW = *maxTime_TPCW_plane1;
				for(int k = 0; k < trackHitTimes_plane1.GetSize(); k++)
				{
					trackHitTimes.push_back(trackHitTimes_plane1[k]);
					trackHitWires.push_back(trackHitWires_plane1[k]);
					trackHitTPCs.push_back(trackHitTPCs_plane1[k]);
					trackHitXvals.push_back(trackHitXvals_plane1[k]);
					trackHitYvals.push_back(trackHitYvals_plane1[k]);
					trackHitZvals.push_back(trackHitZvals_plane1[k]);
					trackHitIsOnTraj.push_back(trackHitIsOnTraj_plane1[k]);
				}
				for(int k = 0; k < trackWf_plane1.GetSize(); k++)
				{
					trackWf.push_back(trackWf_plane1[k]);
				}
			}
			else
			{
				minTime_TPCE = *minTime_TPCE_plane2;
				maxTime_TPCE = *maxTime_TPCE_plane2;
				minTime_TPCW = *minTime_TPCW_plane2;
				maxTime_TPCW = *maxTime_TPCW_plane2;
				for(int k = 0; k < trackHitTimes_plane2.GetSize(); k++)
				{
					trackHitTimes.push_back(trackHitTimes_plane2[k]);
					trackHitWires.push_back(trackHitWires_plane2[k]);
					trackHitTPCs.push_back(trackHitTPCs_plane2[k]);
					trackHitXvals.push_back(trackHitXvals_plane2[k]);
					trackHitYvals.push_back(trackHitYvals_plane2[k]);
					trackHitZvals.push_back(trackHitZvals_plane2[k]);
					trackHitIsOnTraj.push_back(trackHitIsOnTraj_plane2[k]);
				}
				for(int k = 0; k < trackWf_plane2.GetSize(); k++)
				{
					trackWf.push_back(trackWf_plane2[k]);
				}
			}

			int whichTPC = -1;
			if(maxTime_TPCE - minTime_TPCE > maxTime_TPCW - minTime_TPCW)
			{
				whichTPC = 0;
			}
			else
			{
				whichTPC = 1;
			}

			float minT_time = 999999;
			float minT_wire;
			float minT_Xval = 999999;
			float minT_Yval;
			float minT_Zval;
			float maxT_time = -999999;
			float maxT_wire;
			float maxT_Xval = -999999;
			float maxT_Yval;
			float maxT_Zval;
			float minW_low = 999999;
			float maxW_low = -999999;
			float minW_high = 999999;
			float maxW_high = -999999;
			float minW_east = 999999;
			float maxW_east = -999999;
			float minW_west = 999999;
            float maxW_west = -999999;
			float minT_fortimerange = 999999;
			float maxT_fortimerange = -999999;
			for(int k = 0; k < trackHitWires.size(); k++)
			{
				if((trackHitIsOnTraj[k] == true) && (((whichTPC == 0) && (trackHitTPCs[k] == 0)) || ((whichTPC == 1) && (trackHitTPCs[k] == 1))))
				{	
					if((fabs(WF_cathode-fabs(trackHitXvals[k])) > WF_ACdist-maxDistOffset) && (fabs(WF_cathode-fabs(trackHitXvals[k])) < WF_ACdist-minDistOffset) && (fabs(fabs(WF_cathode-fabs(trackHitXvals[k])) - (WF_ACdist-(minDistOffset+maxDistOffset)/2.0)) < fabs(fabs(WF_cathode-fabs(minT_Xval)) - (WF_ACdist-(minDistOffset+maxDistOffset)/2.0))))
					{
						minT_time = trackHitTimes[k];
						minT_wire = trackHitWires[k];
						minT_Xval = trackHitXvals[k];
						minT_Yval = trackHitYvals[k];
						minT_Zval = trackHitZvals[k];
					}
					if((fabs(WF_cathode-fabs(trackHitXvals[k])) > minDistOffset) && (fabs(WF_cathode-fabs(trackHitXvals[k])) < maxDistOffset) && (fabs(fabs(WF_cathode-fabs(trackHitXvals[k])) - (minDistOffset+maxDistOffset)/2.0) < fabs(fabs(WF_cathode-fabs(maxT_Xval)) - (minDistOffset+maxDistOffset)/2.0)))
					{
						maxT_time = trackHitTimes[k];
						maxT_wire = trackHitWires[k];
						maxT_Xval = trackHitXvals[k];
						maxT_Yval = trackHitYvals[k];
						maxT_Zval = trackHitZvals[k];
					}

					//Looks for the minimum wire number that a signal is deposited on
					/*if(((whichTPC == 0) && (trackHitTPCs[k] == 0)) || ((whichTPC == 1) && (trackHitTPCs[k] == 2)))
					{
						if(trackHitWires[k] < minW_low) minW_low = trackHitWires[k];
						if(trackHitWires[k] > maxW_low) maxW_low = trackHitWires[k];
					}
					if(((whichTPC == 0) && (trackHitTPCs[k] == 1)) || ((whichTPC == 1) && (trackHitTPCs[k] == 3)))
					{
						if(trackHitWires[k] < minW_high) minW_high = trackHitWires[k];
						if(trackHitWires[k] > maxW_high) maxW_high = trackHitWires[k];
					}*/

					if(whichTPC == 0 && trackHitTPCs[k] == 0)
                    {
                        if(trackHitWires[k] < minW_east) minW_east = trackHitWires[k];
                        if(trackHitWires[k] > maxW_east) maxW_east = trackHitWires[k];
                    }
					if(whichTPC == 1 && trackHitTPCs[k] == 1)
                    {
                        if(trackHitWires[k] < minW_west) minW_west = trackHitWires[k];
                        if(trackHitWires[k] > maxW_west) maxW_west = trackHitWires[k];
                    }

					if(trackHitTimes[k] < minT_fortimerange)
					{
						minT_fortimerange = trackHitTimes[k];
					}
					if(trackHitTimes[k] > maxT_fortimerange)
					{
						maxT_fortimerange = trackHitTimes[k];
					}
				}
			}

			if((minT_time < 300.0) || (minT_time > maxTicks-300.0)) continue;
			if((maxT_time < 300.0) || (maxT_time > maxTicks-300.0)) continue;

			if((minT_Yval < WF_bottom+10.0) || (minT_Yval > WF_top-10.0)) continue;
			if((maxT_Yval < WF_bottom+10.0) || (maxT_Yval > WF_top-10.0)) continue;

			if((minT_Zval < WF_upstream+10.0) || (minT_Zval > WF_downstream-10.0)) continue;
			if((maxT_Zval < WF_upstream+10.0) || (maxT_Zval > WF_downstream-10.0)) continue;

			if((minT_Zval < 10.0) && (minT_Zval > -10.0)) continue;
			if((maxT_Zval < 10.0) && (maxT_Zval > -10.0)) continue;

			////// STUDIES OF SPATIAL DEPENDENCE //////
			float thewireval = minT_Zval + (minT_Yval - WF_bottom)/sqrt(3);
			//if((thewireval < -700.0) || (thewireval > -620.0)) continue;
			//if((thewireval < -400.0) || (thewireval > -320.0)) continue;
			//if((thewireval < -100.0) || (thewireval > -20.0)) continue;

			float wire_range;
			/*if((fabs(maxW_low-minW_low) > 10000) && (fabs(maxW_high-minW_high) > 10000))
			{
				continue;
			}
			else if(fabs(maxW_low-minW_low) > 10000)
			{
				wire_range = fabs(maxW_high-minW_high);
			}
			else if(fabs(maxW_high-minW_high) > 10000)
			{
				wire_range = fabs(maxW_low-minW_low);
			}
			else
			{
				wire_range = fabs(maxW_low-minW_low) + fabs(maxW_high-minW_high);
			}*/
			if((fabs(maxW_east-minW_east) > 10000) && (fabs(maxW_west-minW_west) > 10000))
            {
                continue;
            }
            else if(fabs(maxW_east-minW_east) > 10000)
            {
                wire_range = fabs(maxW_west-minW_west);
            }
            else if(fabs(maxW_west-minW_west) > 10000)
            {
                wire_range = fabs(maxW_east-minW_east);
            }
            else
            {
                wire_range = fabs(maxW_east-minW_east) + fabs(maxW_west-minW_west);
            }

			float wire_dir_ind1[3] = {0.0, 0.0, 1.0};
			float wire_dir_ind2[3] = {0.0, sqrt(3.0)/2.0, 0.5};
			float wire_dir_col[3] = {0.0, sqrt(3.0)/2.0, -0.5};
			if(whichTPC == 1)
			{
				wire_dir_ind2[2] = -0.5;
				wire_dir_col[2] = 0.5;
			}

			if (runformat == "mc")
			{
				float driftVel = 0.157565;
			}
			
			float track_angle_phi = (180.0/3.14159265)*atan(((maxT_fortimerange-minT_fortimerange)*timeTickSF*driftVel)/(wire_range*wirePitch));
			//if(planenum == 0)
			//{
			//  track_angle_phi = (180.0/3.14159265)*atan((*trackDirX)/((*trackDirY)*wire_dir_ind1[2] - (*trackDirZ)*wire_dir_ind1[1]));
			//}
			//else if(planenum == 1)
			//{
			//  track_angle_phi = (180.0/3.14159265)*atan((*trackDirX)/((*trackDirY)*wire_dir_ind2[2] - (*trackDirZ)*wire_dir_ind2[1]));
			//}
			//else if(planenum == 2)
			//{
			//  track_angle_phi = (180.0/3.14159265)*atan((*trackDirX)/((*trackDirY)*wire_dir_col[2] - (*trackDirZ)*wire_dir_col[1]));
			//}
			//if(track_angle_phi < 0.0)
			//{
			//  track_angle_phi *= -1.0;
			//}
			//cout << track_angle_phi << " " << (180.0/3.14159265)*atan(((maxT_fortimerange-minT_fortimerange)*timeTickSF*driftVel)/(wire_range*wirePitch)) << endl;

			float track_angle_theta;
			if(planenum == 0)
			{
				track_angle_theta = (180.0/3.14159265)*acos((*trackDirY)*wire_dir_ind1[1] + (*trackDirZ)*wire_dir_ind1[2]);
			}
			else if(planenum == 1)
			{
				track_angle_theta = (180.0/3.14159265)*acos((*trackDirY)*wire_dir_ind2[1] + (*trackDirZ)*wire_dir_ind2[2]);
			}
			else if(planenum == 2)
			{
				track_angle_theta = (180.0/3.14159265)*acos((*trackDirY)*wire_dir_col[1] + (*trackDirZ)*wire_dir_col[2]);
			}
			if(track_angle_theta > 90.0)
			{
				track_angle_theta = 180.0-track_angle_theta;
			}

			int anode_minADC = 999999;
			int anode_maxADC = -999999;
			int anode_timeIndex;
			int anode_wireIndex;
			int cathode_minADC = 999999;
			int cathode_maxADC = -999999;
			int cathode_timeIndex;
			int cathode_wireIndex;
			int total_wires = trackWf.size();
			bool anode_flag = false;
			bool cathode_flag = false;
			for(int k = 0; k < total_wires; k++)
			{
				if(((whichTPC == 0) && (trackWf[k].tpc == 0)) || ((whichTPC == 1) && (trackWf[k].tpc == 1)))
				{
					if(trackWf[k].wire == minT_wire)
					{
						anode_flag = true;

						for(int h = 0; h < trackWf[k].adcs.size(); h++)
						{
							if(((planenum < 2) && (trackWf[k].adcs[h] < anode_minADC)) || ((planenum == 2) && (trackWf[k].adcs[h] > anode_maxADC)))
							{
								anode_timeIndex = h;
								anode_wireIndex = k;
								if(planenum < 2)
								{
									anode_minADC = trackWf[k].adcs[h];
								}
								else if(planenum == 2)
								{
									anode_maxADC = trackWf[k].adcs[h];
								}
							}
						}
					}
					else if(trackWf[k].wire == maxT_wire)
					{
						cathode_flag = true;
						for(int h = 0; h < trackWf[k].adcs.size(); h++)
						{
							if(((planenum < 2) && (trackWf[k].adcs[h] < cathode_minADC)) || ((planenum == 2) && (trackWf[k].adcs[h] > cathode_maxADC)))
							{
								cathode_timeIndex = h;
								cathode_wireIndex = k;
								if(planenum < 2)
								{
									cathode_minADC = trackWf[k].adcs[h];
								}
								else if(planenum == 2)
								{
									cathode_maxADC = trackWf[k].adcs[h];
								}
							}
						}
					}
				}
			}

			if((anode_flag == false) || (cathode_flag == false))
			{
				continue;
			}

			int anode_wireIndex_set[11];
			int cathode_wireIndex_set[11];

			int wireSF;
			if(trackWf[anode_wireIndex].channel < trackWf[cathode_wireIndex].channel)
			{
				wireSF = 1;
			}
			else if(trackWf[anode_wireIndex].channel > trackWf[cathode_wireIndex].channel)
			{
				wireSF = -1;
			}
			else
			{
				continue;
			}

			if((anode_wireIndex-11 < 0) || (anode_wireIndex+11 >= total_wires) || (cathode_wireIndex-11 < 0) || (cathode_wireIndex+11 >= total_wires))
			{
				continue;
			}

			anode_wireIndex_set[5] = anode_wireIndex;
			cathode_wireIndex_set[5] = cathode_wireIndex;
			int anode_central_channel = trackWf[anode_wireIndex].channel;
			int cathode_central_channel = trackWf[cathode_wireIndex].channel;	

			int anode_count = 1;
			int cathode_count = 1;
			for(int i = 1; i <= 11; i++)
			{
				if((anode_count < 6) && (trackWf[anode_wireIndex+i].channel == anode_central_channel+anode_count))
				{
					anode_wireIndex_set[5+wireSF*anode_count] = anode_wireIndex+i;
					anode_count++;
				}
				if((cathode_count < 6) && (trackWf[cathode_wireIndex+i].channel == cathode_central_channel+cathode_count))
				{
					cathode_wireIndex_set[5+wireSF*cathode_count] = cathode_wireIndex+i;
					cathode_count++;
				}
			}
			if((anode_count < 6) || (cathode_count < 6))
			{
				continue;
			}

			anode_count = 1;
			cathode_count = 1;
			for(int i = 1; i <= 11; i++)
			{
				if((anode_count < 6) && (trackWf[anode_wireIndex-i].channel == anode_central_channel-anode_count))
				{
					anode_wireIndex_set[5-wireSF*anode_count] = anode_wireIndex-i;
					anode_count++;
				}
				if((cathode_count < 6) && (trackWf[cathode_wireIndex-i].channel == cathode_central_channel-cathode_count))
				{
					cathode_wireIndex_set[5-wireSF*cathode_count] = cathode_wireIndex-i;
					cathode_count++;
				}
			}
			if((anode_count < 6) || (cathode_count < 6))
			{
				continue;
			}

			float temp_anodeWFvals[2*max_offset+1][11][2*timePadding+1];
			float temp_cathodeWFvals[2*max_offset+1][11][2*timePadding+1];

			float anode_tdc0 = trackWf[anode_wireIndex].tdc0;
			float cathode_tdc0 = trackWf[cathode_wireIndex].tdc0;

			float nonlinearity_A = 99999999.0;
			float nonlinearity_C = 99999999.0;

			int anode_offset_index = max_offset;
			int cathode_offset_index = max_offset;

			for(int i = -1*max_offset; i <= max_offset; i++)
			{
				float deviation_A = 0.0;
				float deviation_C = 0.0;

				float num_A = 0.0;
				float num_C = 0.0;

				for(int q = 0; q < 11; q++)
				{
					int neighbor_numADCs = trackWf[anode_wireIndex_set[q]].adcs.size();
					float neighbor_tdc0 = trackWf[anode_wireIndex_set[q]].tdc0;
					float expected_time = timePadding+(((q-5)*wirePitch)/(timeTickSF*driftVel))*tan((3.14159265/180.0)*track_angle_phi);

					for(int k = 0; k < 2*timePadding+1; k++)
					{
						float val = 0.0;
						if((anode_tdc0-neighbor_tdc0+anode_timeIndex-timePadding+k+i >= 0) && (anode_tdc0-neighbor_tdc0+anode_timeIndex-timePadding+k+i < neighbor_numADCs))
						{
							val = trackWf[anode_wireIndex_set[q]].adcs[anode_tdc0-neighbor_tdc0+anode_timeIndex-timePadding+k+i];
						}
						else if(anode_tdc0-neighbor_tdc0+anode_timeIndex-timePadding+k+i < 0)
						{
							val = trackWf[anode_wireIndex_set[q]].adcs[0];
						}
						else if(anode_tdc0-neighbor_tdc0+anode_timeIndex-timePadding+k+i >= neighbor_numADCs)
						{
							val = trackWf[anode_wireIndex_set[q]].adcs[neighbor_numADCs-1];
						}

						temp_anodeWFvals[i+max_offset][q][k] = val;

						if(((planenum == 2) && (val > threshold[planenum])) || ((planenum < 2) && (val < -1.0*threshold[planenum])))
						{
							deviation_A += pow((k-expected_time)*cos((3.14159265/180.0)*track_angle_phi),2);
							num_A += 1.0;
						}
					}
				}

				for(int q = 0; q < 11; q++)
				{
					int neighbor_numADCs = trackWf[cathode_wireIndex_set[q]].adcs.size();
					float neighbor_tdc0 = trackWf[cathode_wireIndex_set[q]].tdc0;
					float expected_time = timePadding+(((q-5)*wirePitch)/(timeTickSF*driftVel))*tan((3.14159265/180.0)*track_angle_phi);

					for(int k = 0; k < 2*timePadding+1; k++)
					{
						float val = 0.0;
						if((cathode_tdc0-neighbor_tdc0+cathode_timeIndex-timePadding+k+i >= 0) && (cathode_tdc0-neighbor_tdc0+cathode_timeIndex-timePadding+k+i < neighbor_numADCs))
						{
							val = trackWf[cathode_wireIndex_set[q]].adcs[cathode_tdc0-neighbor_tdc0+cathode_timeIndex-timePadding+k+i];
						}
						else if(cathode_tdc0-neighbor_tdc0+cathode_timeIndex-timePadding+k+i < 0)
						{
							val = trackWf[cathode_wireIndex_set[q]].adcs[0];
						}
						else if(cathode_tdc0-neighbor_tdc0+cathode_timeIndex-timePadding+k+i >= neighbor_numADCs)
						{
							val = trackWf[cathode_wireIndex_set[q]].adcs[neighbor_numADCs-1];
						}

						temp_cathodeWFvals[i+max_offset][q][k] = val;

						if(((planenum == 2) && (val > threshold[planenum])) || ((planenum < 2) && (val < -1.0*threshold[planenum])))
						{
							deviation_C += pow((k-expected_time)*cos((3.14159265/180.0)*track_angle_phi),2);
							num_C += 1.0;
						}
					}
				}

				if(num_A == 0.0)
				{
					nonlinearity_A = -999.0;
					anode_offset_index = max_offset;
				}
				else if(log10(deviation_A/num_A) < nonlinearity_A)
				{
					nonlinearity_A = log10(deviation_A/num_A);
					anode_offset_index = i+max_offset;
				}

				if(num_C == 0.0)
				{
					nonlinearity_C = -999.0;
					cathode_offset_index = max_offset;
				}
				else if(log10(deviation_C/num_C) < nonlinearity_C)
				{
					nonlinearity_C = log10(deviation_C/num_C);
					cathode_offset_index = i+max_offset;
				}
			}

			//if((nonlinearity_A == -999.0) || (nonlinearity_C == -999.0)) continue;

			//q GOES THROUGH 5 WIRES EITHER SIDE OF THE CLOSEST ANODE WIRE TO THE TRACK (14.5CM AWAY)
			//k IS GOING FROM 0 TO 400 TIME TICKS
			//THIS LOOP COLLAPSES 11 WAVEFORMS INTO A SINGLE 1D WAVEFORM    
			for(int q = 0; q < 11; q++)
			{
				for(int k = 0; k < 2*timePadding+1; k++)
				{
					anode_WFvals[q*(2*timePadding+1)+k] = temp_anodeWFvals[anode_offset_index][q][k];
				}
			}
			for(int q = 0; q < 11; q++)
			{
				for(int k = 0; k < 2*timePadding+1; k++)
				{
					cathode_WFvals[q*(2*timePadding+1)+k] = temp_cathodeWFvals[cathode_offset_index][q][k];
				}
			}

			run = *runNum;
			event = *eventNum;
			cryo = *cryoNum;
			tpc = whichTPC;
			plane = planenum;

			phi = track_angle_phi;
			theta = track_angle_theta;

			anode_T = minT_time;
			anode_W = minT_wire;
			anode_X = minT_Xval;
			anode_Y = minT_Yval;
			anode_Z = minT_Zval;
			anode_NL = nonlinearity_A;

			cathode_T = maxT_time;
			cathode_W = maxT_wire;
			cathode_X = maxT_Xval;
			cathode_Y = maxT_Yval;
			cathode_Z = maxT_Zval;
			cathode_NL = nonlinearity_C;

			trackTree->Fill();
			selected_tracks[planenum]++;
		}
	}

	for(int k = 0; k <= 2; k++)
	{
		cout << "Plane " << k << " Selected Tracks:  " << selected_tracks[k] << endl;
	}

	outfile.cd();
	trackTree->Write();
	outfile.Close();

	timer.Stop();
	cout << "WaveformStudy Runtime:  " << timer.RealTime() << " sec." << endl;

	return 0;
}

void AddFiles(TChain *ch, const char *dir, const char* substr, int minfilenum, int maxfilenum)
{
	TString dirname = (TString) dir;
	//dirname.ReplaceAll("/pnfs/","root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/");

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
				if((counter >= minfilenum) && (counter <= maxfilenum))
				{
					cout << counter << " " << dir << "/" << fname.Data() << endl;

					ch->AddFile(Form("%s/%s/caloskim/TrackCaloSkim",dirname.Data(),fname.Data()));
				}

				counter++;
			}
		}
	}

	return;
}

void AddFilesList(TChain *ch, const char *listname, const char* substr, int minfilenum, int maxfilenum)
{
	ifstream infile(listname);
	string filename;

	int counter = 1;
	while(getline(infile, filename))
	{
		TString fname = (TString) filename;
		if(fname.Contains(substr) && fname.EndsWith(".root"))
		{
			if((counter >= minfilenum) && (counter <= maxfilenum))
			{
				cout << counter << " " << filename << endl;

				ch->AddFile(filename.c_str());
				ch->AddFile(Form("%s/caloskim/TrackCaloSkim",filename.c_str()));
			}

			counter++;
		}
	}

	return;
}

void CheckIfDirectoryExists(std::filesystem::path filepath)
{
	if (!std::filesystem::exists(filepath))
	{
		std::filesystem::create_directories(filepath.parent_path());
	}
	else
	{
		std::cout << "File directory " << filepath << " already exists" << std::endl;
	}
}
