./ProcessWaveformFiles $1 $2 0 1 2
./ProcessWaveformFiles $1 $2 1 1 1
./ProcessWaveformFiles $1 $2 2 1 1
#root -l -b -q ComparePlotsWaveform.C
