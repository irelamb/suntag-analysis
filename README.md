# suntag-analysis
Extraction and analysis of SunTag traces, starting from time-lapse microscopy data and spatial trace annotations (obtained with TrackMate).

## 1) SpotFitting
Bayesian inference of the spot intensity for each trace.

**Input**
    
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

**Input**  

- JSON file contaning the GFP intensity traces for control and perturbed conditions
- Stan model (`.stan`)

**Output** 

`.pkl` file containing the inferred model parameters

Two different models can be used 

- `RB2_Sigma_time.stan`: models where the mature protein intensity (`u`) is a model parameter
- `RB2_Sigma_time_uFixed.stan` : model where the mature protein intensity is a fixed parameter, given as input 
