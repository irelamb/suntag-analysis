# suntag-analysis
Extraction and analysis of SunTag traces, starting from time-lapse microscopy data and spatial trace annotations (obtained with TrackMate).

1) SpotFitting
  Bayesian inference of the spot intensity for each trace.

  Input data:
    
    - time-lapse image
    
    - positions of the spot in each trace (TrackMate output file)
    
  Output:
    
    - csv file for each trace and channel, containing the intensity of each spot.
