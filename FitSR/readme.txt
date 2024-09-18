FitSR

This is a program that takes the output waveform of the ProcessWaveformFile code in WaveformStudy and compares it with a theoretical waveform that is the convolution of the field repsonse (garfield-sbnd-v1.json), electronics response (parametrised) and a smeared gaussian accounting for different angled tracks.

FitSR.py is the python file to run the code, there is an alternative jupyter notebook file in the 'Other' directory, if that's your preference. Along with other field response file from ICARUS.

In order to run this code you'll need to create a python virtual environment, following these instructions:
1. Make python virtual environment
	python -m venv env
2. Activate it
	env/bin/activate
3. Install dependencies using requirements.txt file in 'Other' directory
	pip install -r requirements.txt

To run FitSR.py:
python3 FitSR.py ../WaveformStudy/Results/data/run14608/WFResults_Plane2_X_data.root run14608
	sys.arg(1): Path to input WFResult*.root file
	sys.arg(2): Name of the run, deciding which directory to save the output plots.

If you with to run jupyter notebook then run this:
	python -m ipykernel install --user --name=env
	jupyter notebook
In the notebook, in the navigation menu select Kernel -> Change Kernel -> env to use the virtual environment.

NOTE: before commiting anything make sure you clean the output of the notebooks to avoid gunking up the repository:
nbstripout *.ipynb
