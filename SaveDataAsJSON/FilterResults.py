#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:10:59 2023

@author: Irene

Functions for smoothing data and filtering

"""

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

####

def smooth(y, w, mean=True): 
    
    """ 
    Function for smoothing the signal 
    y = numpy array containing the data to smooth 
    w = the size of the window
    mean = using mean if true, else using median.
    """
    
    # 1. removing nans if present
    mask = ~np.isnan(y)
    y_nonan = y[mask]
    
    if y_nonan.size == 0:
        return y
    
    # 2. applying the smoothing
    
    hw = int((w-1)/2) # half window
    
    y_ext = np.empty((len(y_nonan)+2*hw,)) # longer version of y, were data are repeated at the beginning and the end
    
    y_ext[hw:-hw] = y_nonan
    y_ext[:hw] = y_nonan[0]
    y_ext[-hw:] = y_nonan[-1]
    
    if mean:
        ysmoo = pd.Series(y_ext).rolling(window=w, center=True).mean().to_numpy()
    else:
        ysmoo = pd.Series(y_ext).rolling(window=w, center=True).median().to_numpy()
    
    # 3. getting back the original dimensions (and nans)
    
    res = np.full(y.shape, np.nan)
    res[mask] = ysmoo[hw:-hw]
    
    return res

#%%

# def _where_spot_is_lost(yg, yr):
    
#     mask = yg < 50
#     xg = np.arange(len(yg))[mask]
    
    
#     m = np.nanmean(yr)
#     mask = yr < 0.6 * m
#     xr = np.arange(len(yr))[mask]
    
#     res = np.array([x for x in xg if x in xr])

#     return res


# ###
# def correct_for_spot_loss(yg, yr):
    
#     mask = ~np.isnan(yg)
    
#     yg_nonan = yg[mask]
#     yr_nonan = yr[mask]
    
#     yg_corr = yg_nonan.copy()
            
#     i_lost = _where_spot_is_lost(yg_nonan, yr_nonan)
#     if len(i_lost) > 0:
#         yg_corr[i_lost] = [np.nan]*len(i_lost)
#         yg_corr = pd.Series(yg_corr).interpolate()
    
#     res = np.full(yg.shape, np.nan)
#     res[mask] = yg_corr
    
#     return res

def _butter_lowpass_filter(data, cutoff, fs = 1/20, order = 2):
    """
    
    Function that given the data return the low-pass filter version of it.
    
    Parameters
    ----------
    data : array, float
        intensity over time (it shouldn't contain nans.)
    cutoff : float
        desired cutoff frequency of the filter (Hz).
    fs : float
        sample rate (Hz).
    order : int
        the order of the filter
        

    Returns
    -------
    y : array, float, same shape as data
        smoothed intensity over time

    """
    if np.any(np.isnan(data)):
        raise ValueError('data should not contain NaNs.') 
    
    nyq = 0.5 * fs  # Nyquist Frequency
    
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y


def remove_positive_corr(y, y_ref, cutoff = 1 / 160):
    """
    Function that given two signals, remove from the first the
    fluctuations that are correlated with the fluctuations of the second.

    Parameters
    ----------
    y : array, float
        signal of interest.
    y_ref : array, float
        reference signal, used for correction.
    cutoff : float, optional
        cutoff frequency for the low-pass filter. The default is 1 / 160.

    Returns
    -------
    res : array, float, same shape as y
        corrected signal.

    """
    
    
    res = y.copy()
    
    mask = ~np.isnan(y)
    y_butter = _butter_lowpass_filter(y[mask], cutoff)
    y_ref_butter = _butter_lowpass_filter(y_ref[mask], cutoff)
    
    dy = (y[mask] - y_butter) / y_butter
    dy_ref = (y_ref[mask] - y_ref_butter) / y_ref_butter
    
    dY = np.array([dy, dy_ref])
    
    f = np.mean(dY, axis = 0)
    
    res[mask] = y[mask] - f * y_butter
    
    return res


### Alternative correction ###

def remove_spikes(yg, cutoff = 1 / 120):
    """
    Functions that given a signal yg filters shot noise
    (high frequency spikes/drops from a reference low-pass
     filter signal)

    Parameters
    ----------
    y : array, float
        signal of interest.
    cutoff : float, optional
        cutoff frequency for the low-pass filter. The default is 1 / 120.

    Returns
    -------
    res : array, float, same shape as y
        corrected signal.
    """
    
    
    res = yg.copy()
    
    mask = ~np.isnan(yg)
    ygbutter = _butter_lowpass_filter(yg[mask], cutoff)
    
    dyg = (yg[mask] - ygbutter) / ygbutter
    cond = abs(dyg) > 0.52 # 3 * sigma
    
    ff = dyg.copy()
    ff[~cond] = 0

    res[mask] = yg[mask] - ff * ygbutter
    
    return res


#%%


def compute_starting_value(y, w): # Compute intensity at time 0
    """ 
    Function computing the intial value of y after smoothing.
    """
    
    # As I was doing it previously:
    ysmoo = smooth(y, w=w, mean=False) # this contains nans
    ysmoo_nonan = ysmoo[~np.isnan(ysmoo)] # no nans
    initial_intensity = np.max(ysmoo_nonan[:3]) # first minute # max
    
    return initial_intensity
    
    # I changed it to make it consistent with the way we are computing
    # the run-off time:
    # print("new function")
        
    # ysmoo = smooth(y, w=15, mean=False)

    # mask = ~np.isnan(ysmoo)
    # ysmoo_nonan = ysmoo[mask]
    
    # return ysmoo_nonan[0]


####


def is_lower(y, thr, w):
    """
    Function that returns True if the signal y is always lower than thr,
    after smoothing. False otherwise.
    y : numpy array 
    thr : float
    """
    # moving average
    ysmooth = smooth(y, w)
    # removing nans
    mask = ~np.isnan(ysmooth)
    ysmooth = ysmooth[mask]
    
    return np.all(ysmooth<thr)


####


def checkpoint(df_g, df_r):
    """
    Function that given the data of green (df_g) and
    red (df_r) channels as dataframe returns True if 
    the Track satisfies a selection criterion and
    False otherwise.
    """
    w = 7
    
    sigma_g = df_g["sigma (um)"].to_numpy()
    
    y_r = df_r["I (au)"].to_numpy()
    sigma_r = df_r["sigma (um)"].to_numpy()
    
    return is_lower(y_r, 150, w) and is_lower(sigma_r, 0.23, w) and is_lower(sigma_g, 0.25, w)


####