#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 10:17:52 2021

@author: Irene

small class to store the result of the fit for one track.

v1 : Function added that allows to retrieve a numpy array containing the values of a parameter in all frames.
    
"""

import numpy as np
import pandas as pd


class TrackFit:
    
    def __init__(self, frames=[], fits=[]):
        
        self._dic = {}
        
        for frame, fit in zip(frames, fits):
            
            self._dic[frame] = fit
            
    def __setitem__(self, frame, fit):
        """
        
        Parameters
        ----------
        frame : float
            Frame number.
        fit : Stan fit object
            Result of the fit of one spot.

        Returns
        -------
        None.

        """
        
        self._dic[frame] = fit
        
    def __getitem__(self, frame):
        
        return self._dic[frame]
    
    def frames(self):
        
        return list(self._dic.keys())
    
    def get(self, name, feature):
        
        """
        
        Parameters
        ----------
        name : string
            Name of a parameter (e.g. x0).
        feature : string 
            Feature of the parameter, like mean or sd.

        Returns
        -------
        Numpy array containing the specified parameter feature for all frames.

        """
        
        res = np.empty(len(self.frames()))
        
        for i, frame in enumerate(self.frames()):
            
            fit = self.__getitem__(frame)
            summary_dict=fit.summary()
            df = pd.DataFrame(summary_dict['summary'],
                      columns=summary_dict['summary_colnames'],
                      index=summary_dict['summary_rownames'])
            
            res[i] = df.loc[name, feature]
        
        return res 