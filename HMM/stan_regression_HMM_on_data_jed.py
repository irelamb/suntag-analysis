#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 17:37:52 2023

@author: Irene

regression in Stan to fit the full traces, infer the step size and the number of 
ribosomes per time point

"""

from cmdstanpy import CmdStanModel

from os.path import join

##### IT service #####

from cmdstanpy import cmdstan_path, set_cmdstan_path, show_versions
print(cmdstan_path())

#set_cmdstan_path('/home/ilambert/.cmdstan/cmdstan-2.34.1')

print(cmdstan_path())
print(show_versions())

#####

#import stan
import numpy as np
import json
from os.path import join

#import nest_asyncio
#nest_asyncio.apply()


import pickle

from sys import argv

homedir = __file__ # where the compiled models (.stan) will be saved

## to be modified by the users
workingdir = "path_to_results" # where the inference results will be saved
inputdir = "path_to_JSON_input_file"
#####


from time import time

#%%

if __name__ == '__main__':
    
    HMM_model_name = argv[1]
    data_file = argv[2]
    b0 = float(argv[3])
    sigma0 = float(argv[4])
    L = int(argv[5]) # total length of the reporter
    max_step_size = int(argv[6])
    N_max = int(argv[7])
    n_samples = int(argv[8])
    #u = int(argv[7])

    print("Input Paramters")
    print("model name", HMM_model_name)
    print("file", data_file)
    print("L", L)
    print("max_step_size", max_step_size)
    print('\n')
    
    # ---- parameters ----- #
    
    nframes = 271 # number of frames
    dt = 20 # sec
    t = np.arange(nframes) * 20 # time is seconds

    S = 599 # length of the SunTag only (without the last 6AA linking it to the gene )

    #N_max = 30 # increase
    particle_size = 10

    print("nframes", nframes)
    print("dt", dt)
    print("L", L)
    print("S", S)
    print("N_max", N_max)
    print("particle size", particle_size)

    print('\n')
    print("n_samples", n_samples)

    #-------- Load data ----------#
    
    u = 14.0 # from calibration experiments


    json_file = join(inputdir, data_file)
    with open(json_file) as f:
        data = json.load(f)
    
    data['b0'] = b0
    data['sigma0'] = sigma0

    data['N_max'] = N_max
    data['max_step_size'] = max_step_size
    data['l'] = particle_size
    data['S'] = S
    data['L'] = L
    data['u'] = u
#    data['sigma'] = sigma


    #----- Compile the model and fit ------#
    
    start = time()
    compiled_model = CmdStanModel(stan_file = join(homedir, '{}.stan'.format(HMM_model_name)))#, force_compile = True)
    FIT = compiled_model.sample(data = data, chains = 2, parallel_chains = 2)
    end = time()
    print("model took ", end-start)

    #----- Save fit ------#
    
    # save the fit
    file = join(workingdir, "{}_{}_N{}_b0{}_s0{}_u{}_step{}.pkl".format(HMM_model_name, data_file[: -5], N_max, b0, sigma0, u,  max_step_size))
    with open(file, "wb") as f:
        pickle.dump({'model' : compiled_model, 'fit' : FIT}, f, protocol=-1)
    print("model saved as ", file)
