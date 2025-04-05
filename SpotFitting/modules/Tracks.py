#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 17:43:09 2021

@author: Irene

a small class to store tracking data from trackmate
"""

import numpy as np
import pandas as pd

from os.path import join


class Tracks:
    
    def __init__(self, path_to_csv):
        """ class with two attributes: tracks and tracks_noNan """
        
        columns = ['FRAME', 'TRACK_ID', 'POSITION_X', 'POSITION_Y', 'ID']
        table = pd.read_csv(path_to_csv, index_col=0, usecols=columns)
        
        # private attributes
        self._tracks = self._dataframe_to_dictionary(table)
        self._replace_nans()
        
        # public attributes
        self.IDs = list(self._tracks.keys())
        
        
    def __getitem__(self, trackID):
        
        return self._tracks[trackID]
        
            
    def _dataframe_to_dictionary(self, table):
        
        """ Convert the dataframe to a dictionary with track ID as keys 
        and (t,x,y) dataframe as values """
        
        dic = {}
        
        trackIDs = np.unique(table["TRACK_ID"])
        
        for ID in trackIDs:
            bool_select = table['TRACK_ID']==ID
            T = table[bool_select]
            T = T.set_index('FRAME', drop=True)
            T = T.sort_index()
            dic[ID] = T.drop('TRACK_ID', axis=1)
        
        self._add_Nans(dic)
        
        return dic
    
    
    def _add_Nans(self, dic):
        """ Add Nans for frames with missing detection. """
        
        for key, df in dic.items():
            
            min_frame = df.index[0]
            max_frame = df.index[-1]
    
            new_index = np.arange(min_frame, max_frame+1)
            
            df = df.reindex(new_index) # this will insert NaNs
            
            dic[key] = df
            
            
    def _replace_nans(self):
        
        """ replace Nans with valid values. """
        
        print("Replacing missing detections with the last valid position...")
        
        for ID, track in self._tracks.items():
            
            missing = np.isnan(track["POSITION_X"])
            missing_frames = track.index[missing]

            while len(missing_frames):
                
                track.loc[missing_frames, "POSITION_X"] = \
                track.loc[missing_frames - 1, "POSITION_X"].to_numpy()
                
                track.loc[missing_frames, "POSITION_Y"] = \
                track.loc[missing_frames - 1, "POSITION_Y"].to_numpy()
                
                missing = np.isnan(track["POSITION_X"])
                missing_frames = track.index[missing]
        
        print("Done!")          
