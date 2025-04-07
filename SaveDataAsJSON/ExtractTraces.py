#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed April  3 15:58:09 2023

@author: Irene

Functions to load data from files.

v2) in function get_values() append to lists and then convert to ndarrays, instead
    of initializing an array of fixed size. In this way, the shape of the result
    will correspond to the number of traces effectively selected.
    
v3) add stage argument to get_files_list; 
    get_values returns also the array of track ids
    
v4) get_values also returns the array of stages

v5)

v6) if len(y_g) < 15:
        continue
    Added to exclude tracks that are too short (shorter than 5 minutes)
    because I had some issues with a newer trackmate version
    and too short traces were included

"""

import numpy as np
import pandas as pd
from os.path import join
from glob import glob
import re

from FilterResults import checkpoint, smooth

####


def get_files_list(date, sample, workingdir, model_name, spec="", stage="s*", chan=0):
    
    """
    Description:
        Given some parameters defining a data set, it returns the list of files 
        containing the fit.
    
    Arguments:
        date : (str) experiment date, e.g. '20230104'
        sample : (str) sample name, e.g. 'COL-H_01'
        workingdir : (str) path to fit data
        model_name : (str) name of the model used for the fit
        spec : (str) ''(default) or '_ff' for flat-field corrected images
        stage : (str) "s1", "s2", "s3", "s4" or s* (default) if all stages
        chan : (int) 0 (default) or 1
    
    Returns:
        list of files
    """

    # 1. Get some text strings for file names
    
    sample_nonum = sample.split("_")[0]
    sample_num = sample.split("_")[-1]
    
    filename_noext = sample_nonum + "_" + date + "_" + sample_num + "_" + stage
        
    # 2. Find all fit files
    
    files = glob( join(workingdir, date, sample, model_name, 'log*-track*-chan{}'.format(chan), '{}_track*_ch{}_{}{}.csv'.format(filename_noext, chan, model_name, spec)) )   
    print( join(workingdir, date, sample, model_name, 'log*-track*-chan{}'.format(chan), '{}_track*_ch{}_{}{}.csv'.format(filename_noext, chan, model_name, spec)) )
   
    return files


####

def get_trackIDs(files):
    
    """
    return the track ID from a track file.
    """
    
    res = np.array([int(re.search("track(\d+)", f)[1]) for f in files])
    
    return res

####

def get_stages(files):
    
    """
    return the stage number from a track file.
    """
    
    res = np.array([int(re.search("_s(\d)", f)[1]) for f in files])
    
    return res

####

def get_dates(files):
    
    """
    return the date (as int) from a track file.
    """
    
    res = np.array([int(re.search("_(\d{8})_", f)[1]) for f in files])
    
    return res 


###

def get_unique_IDs(files):
    
    """
    return an ID for each track that is unique in the same acquisition.
    The ID is composed by the stage number followed by the track number.
    """
    
    unique_IDs = []
    
    for f in files:
        
        ID = re.search('track(\d+)', f).group(1)
        stage = re.search('_s(\d)', f).group(1)
        
        unique_IDs.append(int(stage+ID))
        
    return np.array(unique_IDs)


# def get_values(files, nframes, column='I (au)'):
    
#     """
#     Function returning a ndarray of dimensions (len(files), nframes)
#     containing the values of column (default intensity).
#     """
    
#     res = np.full( (len(files), nframes), np.nan )
    
#     for n, file in enumerate(files):
        
#         try:
#             df = pd.read_csv(file, index_col=0)
#         except FileNotFoundError as e:
#             print(e)
#             continue
        
#         y = df[column].to_numpy()
#         # time
#         t = df["time (s)"].to_numpy()
#         dt = t[1]-t[0]
#         index = (t//dt).astype(int)
        
#         #
        
#         res[n][index] = y
    
#     return res


####


def get_values(files_g, nframes, column_name="I (au)"):
    
    """
    Description:
        Given a list of .csv files containing the result of the fit in the green
        channel, it retrieve the values of column_name both for red and green 
        channel. 
    Parameters:
        files_g : (list) list of .csv files with the green channel fit.
        nframes : (int) maximum number of frames per track.
        column_name : (str) name of the desired measurement.
    Returns:
        traces_g : nd-array containing the green traces, each of length nframes.
        traces_r : nd-array containing the red traces, each of length nframes.
        mask : a mask for the traces with a red signal that satisfies some filtering
        criteria.
    """
    
    traces_g = []
    traces_r = []
    new_files = [] # new file list corresponding to the traces
    # excluding files that do not have a correspondingred trace
    
    # to do the filtering
    I_g = []
    I_r = [] # intensity red channel
    W_r = [] # size of the spot in the red channel
    
    for n, file in enumerate(files_g):
        
        try:
            # green channel
            df_g = pd.read_csv(file, index_col=0)
            
            # red channel
            file_r = re.sub("ch0", "ch1", file)
            file_r = re.sub("chan0", "chan1", file_r)
            df_r = pd.read_csv(file_r, index_col=0)
            
            
        except FileNotFoundError as e:
            print(e)
            continue
        
        
        
        y_g = df_g[column_name].to_numpy()
        y_r = df_r[column_name].to_numpy()
        
        if len(y_g) < 15: # traces longer than 5 min.
            continue
        
        # time
        t = df_g["time (s)"].to_numpy()
        dt = t[1]-t[0]
        index = (t//dt).astype(int)
        
        #
        temp = np.full( nframes, np.nan )
        temp[index] = y_g
        traces_g.append(temp)
            
        temp = np.full( nframes, np.nan )
        temp[index] = y_r
        traces_r.append(temp)
        
        # for the filtering
        i_g = df_g['I (au)'].to_numpy()
        i_r = df_r['I (au)'].to_numpy()
        w_r = df_r['sigma (um)'].to_numpy()
        
        temp = np.full( nframes, np.nan )
        temp[index] = i_g
        I_g.append(temp)
        
        w=7
        temp = np.full( nframes, np.nan )
        temp[index] = i_r
        I_r.append(smooth(temp,w))
        
        temp = np.full( nframes, np.nan )
        temp[index] = w_r
        W_r.append(smooth(temp,w))
        
        new_files.append(file)
    
    new_files = np.array(new_files)
    traces_g = np.array(traces_g, dtype = float)
    traces_r = np.array(traces_r, dtype = float)
    I_g = np.array(I_g, dtype = float)
    I_r = np.array(I_r, dtype = float)
    W_r = np.array(W_r, dtype = float)
    
    #print(len(I_g), len(I_r), len(W_r))
    
    
    median = np.nanmedian(I_r)
    mask =  np.array([np.all(y[~np.isnan(y)] < 2 * median) for y in I_r], dtype = bool) & \
            np.array([np.nanmean(y) < 100 for y in I_r], dtype = bool) # this depends on the settings of the experiment

    
    mask = mask & \
           np.array([np.all(y[~np.isnan(y)] < 0.23) for y in W_r], dtype = bool) # spot width
           
    mask = mask & \
           np.array([ np.all(np.abs(y[~np.isnan(y)][1:] - y[~np.isnan(y)][:-1]) < 100) for y in I_g]) # very big drops or increases
    
    return new_files, traces_g, traces_r, mask
