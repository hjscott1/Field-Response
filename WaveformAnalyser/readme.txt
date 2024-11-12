WaveformStudy

This code package uses the output histroot Reco2 files from caloskim ntuples with RawDigits to select anode-cathode crossing cosmic muon tracks, 2-5cm from the anode, and reconstruct the average output signal on TPC APA readout wires. The code is separated into multiple steps:

1. WaveformStudy: Provides the initial production of waveform information, now stored as an ntuple.
	Example: ./WaveformStudy InputFiles/WireBias1.list 1 10 data test Reco2
	 - arg[1]: filepath - Directory with calibration ntuples in
	 - arg[2]: filenums - Number files to run through in the filepath
	 - arg[3]: runformat - mc or data
	 - arg[4]: outputfilename - Name of folder to store the results
	 - arg[5]: filetext - (optional, default “root”), Runs all files in arg[1]’s directory with this pattern
	Output: wfstudy_data.root or wfstudy_mc.root

2. ProcessWaveformFiles: Takes the individual waveforms from the ntuple and combines them together.
	Example: ./ProcessWaveformFiles ./Results/data/run14608/wfstudy_data.root X 2 1 1
	 - arg[1]: filepath - Directory that WaveformStudy output in
	 - arg[2]: runtype - Chooses tracks in specific cryostats and TPCs. Can be X (both TPCs, -1), E (east TPC, 0), W (west TPC, 1).
	 - arg[3]: planenum - Chooses which wire plane to select signals from. 1st induction: 0, 2nd induction: 1, collection: 2.
	Output: WFresults_Plane2_X_data.root

3a. MakePlotsWaveform.C: Compares mc and data waveforms, smearing the mc waveforms to account for different angled tracks (will need some edits to automate it).
	Example: root -l -b -q MakePlotsWaveform.C

3b. ComparePlotsWaveform.C - Directly compares any two waveforms
	Example: root -l -b -q ComparePlotsWaveform.C

3c. Compare3PlotsWaveform.C - Directly compares any three waveforms
	Example: root -l -b -q Compare3PlotsWaveform.C

3d. CompareTPCsWaveform.C - Directly compares the waveforms from the east, west and resultant TPCs.
	Example: root -l -b -q CompareTPCsWaveform.C
