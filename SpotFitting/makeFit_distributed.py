#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 17:29:07 2021

@author: Irene

1) using physical units instead of pixel units for the fit
2) add line to check if track fit already exists.
3) arguments as DATE, REP, CHAN and STAGE are passed as command line arguments, in this order. If no STAGE is passed the variable is assigned to empty string.
4) params_tocsv function added, to save some parameters to a csv file.
5) read parameters directly as inputs.
6) (2023-02-07) adding to the output csv file all the parameters
7) (2024-09-04) v5: regarding the prior on i0 (background intensity): because the background is highly variable,
    especially in the red channel, the idea is to set the mean of i0 prior to the median values of the outer
    pixels in a spot frame.
8) use the outer pixels to compute both the mean and the standard deviation of the i0 prior.
"""

#from params import workingdir, homedir, model_name, \
#    cmap, cyto_avg, cyto_std

homedir = "/home/ilambert/SunTag"
workingdir = "/scratch/ilambert/SunTag/"

channels = ["sdcGFPquad", "sdcCy5quad"] # channels names
cmap = ["Greens_r", "Reds_r"] # colors for plotting
       
import sys
sys.path.append(homedir) # necessary to import the modules

# Import Image and Track modules
from modules.MicroscopyImage import MicroscopyImage
from modules.Tracks import Tracks
from modules.TrackFit import TrackFit

# Import Python packages
import numpy as np
from math import ceil
import pystan
import time

import re

# for plotting
#import matplotlib.pyplot as plt

import pickle
import csv
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')

from os.path import join
import os
from sys import argv

import multiprocessing as mp

import modules.models as models

def is_ometif(filename):
    l = filename.split('.')
    return len(l)==3

# -------------------------------------------------------------------------- #


# Spot-fitting function
def fit_spot(model, intensities, xpos, ypos, x0, y0, bg0, bg0_std):
    
    """
    Function calling Stan HMC sampling method.

    - intensities = 2D array containing the intensity values at a red spot position
    - xpos, ypos = x,y pixel coordinates in physical units
    - x0, y0 = position in physical unit of the red spot center
    - bg0 = value used as mean intensity value for the background prior (i0)
    - bg0_std = standard deviation of the background prior on i0

    Returns a Stan fit object.
    """
    
    data = {
            'n' : len(xpos),
            'x' : xpos,
            'y' : ypos,
            'intensity' : intensities,
            'x_red' : x0,
            'y_red' : y0,
            'cyto_avg' : bg0,
            'cyto_std' : bg0_std # this is used only in gaussian models
            }

    fit = model.sampling(data=data, iter=1000, chains=1, warmup=500, thin=1, seed=101,
                             verbose=False)
    
    return fit

# -------------------------------------------------------------------------- #

# Log the fits for a track
def rawLog_fits(fit, trackID, frame, n, filename):
    # TESTING: time how long it takes to log the fit
    # timer=time.time()
    
    path2fit = join(workingdir, "Results", DATE, REP, MODEL_NAME, "log-{}track{}-chan{}".format(STAGE, trackID, CHAN))

    csvname="log-{}_track{}_ch{}_{}{}.csv".format(filename, trackID, CHAN, MODEL_NAME, suffix)

    # 2. Get fit summary
    summary_dict=fit.summary()
    df = pd.DataFrame(summary_dict['summary'],
                      columns=summary_dict['summary_colnames'],
                      index=summary_dict['summary_rownames'])

    #3. Log the fit into csv
    if n == 0: # overwrite previous log, if it exists
        with open(join(path2fit, csvname), 'w+') as log:
            log_w = csv.writer(log, delimiter=',')
            log_w.writerow(["dt = {}, {}, chan = {}, {}".format(dt, FILENAME, CHAN, MODEL_NAME)])
            log_w.writerow(["Frame {}".format(frame)])
    else: #append following records to the log
        with open(join(path2fit, csvname), 'a') as log:
            log_w=csv.writer(log, delimiter=',')
            log_w.writerow(["Frame {}".format(frame)])
    df.to_csv(join(path2fit, csvname), mode='a')

    # TESTING: time how long it takes to log the fit
    # print("Logging the fit took {}s".format(time.time() - timer))

    return

# -------------------------------------------------------------------------- #

def params_tocsv(fit, trackID, frames, filename):
    
    path2fit = join(workingdir, "Results", DATE, REP, MODEL_NAME, "log-{}track{}-chan{}".format(STAGE, trackID, CHAN))
    csvname="{}_track{}_ch{}_{}{}.csv".format(FILENAME_NOEXT, trackID, CHAN, MODEL_NAME, suffix)
    
    x0 = fit.get(name="x0", feature="mean")
    y0 = fit.get(name="y0", feature="mean")
    I = fit.get(name="I", feature="mean")
    sd_I = fit.get(name="I", feature="sd")
    sigma = fit.get(name="sigma", feature="mean")
    sd_sigma = fit.get(name="sigma", feature="sd")
    
    i0 = fit.get(name="i0", feature="mean")
    
    try:
        a = fit.get(name="a", feature="mean")
        b = fit.get(name="b", feature="mean")
        sigma_ln = fit.get(name="sigma_ln", feature="mean") # in the lognormal model
    except Exception:
        a = np.nan
        b = np.nan
        sigma_ln = np.nan
    
    lambda_ = fit.get(name="lambda", feature="mean")
    tau = fit.get(name="tau", feature="mean")
    
    df = pd.DataFrame({'time (s)' : frames*dt, 'x0 (um)' : x0, 'y0 (um)' : y0, 'I (au)' : I, 'sd_I (au)' : sd_I, 'sigma (um)' : sigma, 'sd_sigma (um)' : sd_sigma, 'i0 (au)' : i0, 'a' : a, 'b' : b, 'lambda' : lambda_, 'tau' : tau, 'sigma_ln' : sigma_ln }, index=frames)
    df.to_csv(join(path2fit, csvname), index_label='frame')
    
    return 
    

# -------------------------------------------------------------------------- #


#Get plots the fits for a track
def plot_fits(fit, trackID, frame, image, intensities, xpos, ypos, filename):
    # TESTING: time how long it takes to make and save a plot
    # timer=time.time()

    # 1. Make csv file name and path

    path2plot = join(workingdir, "Plots", DATE, REP, MODEL_NAME, "log-{}track{}-chan{}".format(STAGE, trackID, CHAN))

    if not os.path.exists(path2plot):
        os.makedirs(path2plot)

    plotname = "plot-{}_track{}_ch{}_{}_frame{}{}.jpg".format(filename, trackID, CHAN, \
                                                                    MODEL_NAME, frame, suffix)

    #2. Create a figure, do colour normalisation

    p10_g, p99_g = np.quantile(image[:, 0], (0.1, 0.999999))
    norm = matplotlib.colors.Normalize(vmin=p10_g, vmax=p99_g)

    fig = plt.figure(figsize=(16, 12))

    #3. Show the image
    ax_im = fig.add_subplot(1, 2, 1)
    ax_im.imshow(intensities, cmap=cmap[CHAN], norm=norm)
    ax_im.set_title("Track {}, frame {}, time {:.2f}min image".format(trackID, frame, frame*dt/60))

    #4. Plot the fit
    #get fit summary
    summary_dict = fit.summary()
    df = pd.DataFrame(summary_dict['summary'],
                      columns=summary_dict['summary_colnames'],
                      index=summary_dict['summary_rownames']) 

    xv, yv = np.meshgrid(xpos, ypos, sparse=False, indexing='ij')

    ax_plot = fig.add_subplot(1, 2, 2)
    surface = get_fit_surface(xv, yv, df['mean'])
    norm = matplotlib.colors.Normalize(vmin=p10_g, vmax=p99_g)
                                                                    
    ax_plot.imshow(surface, cmap=cmap[CHAN], norm=norm)

    ax_plot.set_title("Track {}, frame {}, time {:.2f}min image".format(trackID, frame, frame*dt/60))

    fig.tight_layout()
    #Save the figure
    fig.savefig(join(path2plot,plotname))

    # TESTING: time how long it takes to make and save a plot
    # print("Plotting the fit and saving it took {}s".format(time.time()-timer))

    return

# -------------------------------------------------------------------------- #

#Auxiliary function plotting a single fit # THIS SHOULD BE CHANGED!
def get_fit_surface(x, y, p):
    try:
        s = p['a']*x + p['b']*y + p['i0'] + p['I']/(2*np.pi*p['sigma']**2) * np.exp( - (x-p['x0'])**2 / (2*p['sigma']**2) - (y-p['y0'])**2 / (2*p['sigma']**2) )
    except:
        s = p['i0'] + p['I']/(2*np.pi*p['sigma']**2) * np.exp( - (x-p['x0'])**2 / (2*p['sigma']**2) - (y-p['y0'])**2 / (2*p['sigma']**2) )
    return s

# -------------------------------------------------------------------------- #


def fitTrack(trackID):

    path2fit = join(workingdir, "Results", DATE, REP, MODEL_NAME, "log-{}track{}-chan{}".format(STAGE, trackID, CHAN))

    if not os.path.exists(path2fit):
        os.makedirs(path2fit)
        
    # check if the track was already analysed 
            
    fitFile = '{}_track{}_ch{}_{}{}.pkl'.format(FILENAME_NOEXT, trackID, CHAN, MODEL_NAME, suffix)
        
    if os.path.exists(join(path2fit,fitFile).strip()):
        print("Fit of track {} already exists. Skipping...".format(trackID))
        return
        
    # Selecting the track

    track = tracks[trackID]
    frames = track.index
    
    track_too_close = False # a flag in case a spot in a track is too close 
                            # to the image border
        
    # Initializing an object to store the fit results
        
    trackRes = TrackFit()
    
    # Time how much the fits took
    start = time.time()
    
    for n, frame in enumerate(frames):
            
        x = track.loc[frame, "POSITION_X"]
        y = track.loc[frame, "POSITION_Y"]
            
        i = int(x/dx) 
        j = int(y/dy) 
            
        sub = image[frame, CHAN].transpose() # transposed to have x first and y second
            
        intensities = sub[ i-R : i+R+1, j-R : j+R+1]
    
        xpos = np.arange(i-R, i+R+1) * dx
        ypos = np.arange(j-R, j+R+1) * dy

        outer_pixels = np.concatenate((intensities[0, :], intensities[-1, :], intensities[1:-1, 0], intensities[1:-1, 1])) # pixels bordering the square


        # used to estimate the average background
        bg0 = np.mean(outer_pixels)
        bg0_std = np.std(outer_pixels)
        print("trackID", trackID, "\n")
        print("frame", n, "\n")
        print("average background estimate: ", bg0, "\n")
        print("standard deviation: ", bg0_std, "\n")
        #print("outer pixels", outer_pixels, "\n")
        # one could check if there are "many" zeros, that could arise if the spot is
        # close to the border of the regressed image

        # Running the fit 

        try:
            trackRes[frame] = fit_spot(sm, intensities, xpos, ypos, x, y, bg0, bg0_std)
        except RuntimeError as e:
            print(e)
            print("Spot too close to border:\n")
            print("position in pixels: ", i, j, "\n")
            print("radius: ", R, "\n")
            print("image shape: ", sub.shape, "\n")
            print("Skipping track ", trackID)
            track_too_close = True
            break
        
        # Log fit summary and plot fits for each frame
        rawLog_fits(trackRes[frame], trackID, frame, n, FILENAME_NOEXT)
        plot_fits(trackRes[frame], trackID, frame, image, intensities, xpos, ypos, FILENAME_NOEXT)
        
        
    if track_too_close: 
        track_too_close=False
        return    
    
    if not os.path.exists(path2fit):
        os.makedirs(path2fit)
            
    with open(join(path2fit, fitFile), 'wb') as f:
        pickle.dump(trackRes, f)
        
    
    # Save parameter values to csv file
    params_tocsv(trackRes, trackID, frames, FILENAME_NOEXT)
        
    end = time.time()
    
    delta = end-start
    print("track {0} took {1:.3f} seconds".format(trackID, delta))
    
# -------------------------------------------------------------------------- #
# -------------------------------------------------------------------------- #

# Fit tracks
if __name__ == '__main__':
    
    # Parameters
    
    FILENAME = os.path.basename(argv[1]) # path
    CHAN = int(argv[2])
    #cyto_avg = float(argv[3]) # cytoplasm average intensity # v5: computed for each spot
    #cyto_std = float(argv[3]) # cytoplasm standard deviation
    MODEL_NAME = argv[3]
    DX = float(argv[4]) # pixel size (micron)
    DY = float(argv[5]) # pixel size (micron)
    DT = float(argv[6]) # time step (s)
    
    print(FILENAME, CHAN, MODEL_NAME, DX, DY, DT)
    
    
    # 1. Check whether the file is a OME.TIF
    
    if is_ometif(FILENAME):
        FILENAME_NOEXT = FILENAME[:-8]
        suffix=""
        
    else:
        FILENAME_NOEXT = FILENAME[:-4] 
        
        FILENAME_NOEXT = FILENAME_NOEXT[:-3] # take out the _ff suffix
        suffix="_ff"
    
    # 2. Extract replicate, date and stage (if present) from file name
    
    filename_split = FILENAME_NOEXT.split('_')
    REP_NONUM = filename_split[0]
    DATE = filename_split[1]
    REP_NUM = filename_split[2]
    
    match_stage = re.search('s\d', FILENAME_NOEXT)
    if match_stage:
        STAGE = match_stage.group(0) + '-'
    else:
        STAGE = ""
    
    REP = REP_NONUM+'_'+REP_NUM
   
    
    # 3. Import the specified model from module 'models' and store it as a variable 'model'
    model = getattr(models, MODEL_NAME)

    
    # 4. Load the image
    path2image = join(workingdir, "Data/Images/{}/{}".format(DATE,REP), FILENAME)
    
    image = MicroscopyImage(path2image, DX, DY, DT, channels) 
    print(image.properties)
    
    # 5. Load the csv file containing the tracking results

    path2csv = join(workingdir, "Data/Tracks/{}/{}".format(DATE, REP),
                "{}_Spots.csv".format(FILENAME_NOEXT))
    print(path2csv)
    
    tracks = Tracks(path2csv) # loading tracks
    
    # 6. Define some constants (coversion factors, radius of the square,
    #    channel to fit)
    
    r = 1.  # radius of the square around the spot
            # defining the subset of pixels to be fitted
            # (radius in physical units (um))

    dx = image.properties.dx
    dy = image.properties.dy 
    dt = image.properties.dt
    R = int(ceil(r/dx)) # radius in pixels
        
    # 5. Make a model and pickle it
    
    path2model = join(workingdir, 'Models')
    modelFile = '{}.pkl'.format(MODEL_NAME)
    
    if not os.path.exists( join(path2model, modelFile) ):
        
        # 5a. Compiling and saving the model
    
        sm = pystan.StanModel(model_code = model)
    
        with open(join(path2model, modelFile), 'wb') as f:
            pickle.dump(sm, f)
    
    else: 
        # 5b. Loading the model
        
        with open(join(path2model, modelFile), "rb") as fp:
            sm = pickle.load(fp)
     
    # distributing
    output = mp.Queue() # define an output queue
    processes = [ mp.Process(target = fitTrack, args = (trackID,)) for trackID in tracks.IDs ]

    for p in processes:
        p.start()

    for p in processes:
        p.join()
        
        
