#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 17:26:39 2021

@author: Irene

a small class to store one fluorescence microscopy image (4D)
as numpy array together with the metadata.

2023-01-25

The input image can be an OME.TIF file or a TIF file.
In case of an OME.TIF file, we try to read metadata from the image.
If the image is TIF or if the OME.TIF does not contain metadata, we use
the metadata specified in the params file.

The number of channels and the channels names are always read from the 
params file because in my ome.tif data they are wrong.

2023-01-26
instead of reading a parameter file the metadata are passed as arguments of the
function __init__

"""

# to import image
from skimage import io

# to read metadata
from PIL import UnidentifiedImageError
from PIL import Image as Im
import xml.etree.ElementTree as ET

from collections import namedtuple

from os.path import join, basename

from time import time

class MicroscopyImage:
    """ Class to store the data collected in a fluorescence
    microscopy experiment (intensity and metadata)
    """
    
    Metadata = namedtuple('Metadata', ['dx', 'dy',
                          'nframes', 'dt', 'nchan',
                          'chan_names']) # namedtuple Metadata class
                                         # (it is unmutable)
                                         # maybe a mutable object would be better
    
    def __init__(self, path2image, dx, dy, dt, channels):
        
        print("Loading image...")
        start = time()
        
        self._image = io.imread(path2image)
        nframes = len(self._image) # number of frames
        
        if self._is_ometif(path2image):
            
            try:
            
                self.properties = self._read_metadata(path2image, channels) 
            
            except UnidentifiedImageError:
                
                print("OME.TIF image does not contain metadata.")
                print("Using parameters file instead.")
                
                self.properties = self.Metadata(dx, dy, nframes, dt, len(channels), channels) 
        
        else: # the image does not contain metadata, use metadata in params file.
        
            self.properties = self.Metadata(dx, dy, nframes, dt, len(channels), channels) 
        
        
        end = time()
        print("Done! (in {}s)".format(end-start))
        
    def __getitem__(self, select):
        
        return self._image[select]
    
    def _is_ometif(self, path2image):
        
        image_name = basename(path2image)
        l = image_name.split('.')
        return len(l)==3
        
        
    def _read_metadata(self, path_to_ometif, channels):
        """ Read metadata from image file.
        at the moment the image can be 4d,
        with x,y,t and channel as dimensions. 
        """
        
        file = Im.open(path_to_ometif) # lazy operation
        
        description = file.tag[270][0]
        
        root = ET.fromstring(description)

        img = root.find("{http://www.openmicroscopy.org/Schemas/OME/2016-06}Image")
        pixels = img.find("{http://www.openmicroscopy.org/Schemas/OME/2016-06}Pixels")

        
        dx = float(pixels.attrib['PhysicalSizeX'])
        dy = float(pixels.attrib['PhysicalSizeY'])
        nframes = int(pixels.attrib['SizeT'])
        dt = float(pixels.attrib['TimeIncrement'])
        
        #channels = img.attrib['Name'].split('/') # this is usually wrong
        
        return self.Metadata(dx, dy, nframes, dt, len(channels), channels)
    
               
if __name__ == "__main__":
    
    PATH = "/Volumes/Naef-Lab/Irene/SunTag_svfas5"
    
    path2ome = join(PATH, "data/corrected/20230113/COL-H-GC7-1uM-24h_01", "COL-H-GC7-1uM-24h_20230113_01_s1_ff.tif")
    
    image = MicroscopyImage(path2ome)
    
    print(image.properties)
    
    io.imshow(image[10,1,:,:]) # dimensions : time, channel, y, x 
        
