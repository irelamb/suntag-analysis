# suntag-analysis
Extraction and analysis of SunTag traces, starting from time-lapse microscopy data and spatial trace annotations (obtained with TrackMate).

## 1) SpotFitting
Bayesian inference of the spot intensity for each trace.

**Input data**
    
- time-lapse image
    
- positions of the spot in each trace (TrackMate output file)
    
**Output**

- csv file for each trace and channel, containing the intensity of each spot.

To run on a server in a distributed way (multiple acquisitions):

`bash 01_submit_job.sh`

## 2) SaveDataAsJSON
Jupyter notebook that merges the files containing the trace intensities into one unique JSON file that will be the input of Bayesian inference.

## 3) HMM
Bayesian inference with Hidden Markov Model (HMM)

**Input data**  JSON file contaning the GFP intensity traces

**Output data** .pkl file containing the inferred model parameters
