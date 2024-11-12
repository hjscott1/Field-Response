./ProcessWaveformFiles $1 $2 0 
./ProcessWaveformFiles $1 $2 1
./ProcessWaveformFiles $1 $2 2
#root -l -b -q ComparePlotsWaveform.C
