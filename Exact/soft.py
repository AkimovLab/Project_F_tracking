import sys
import cmath
import math
import os
import h5py
import matplotlib.pyplot as plt   # plots
import numpy as np
import time
import warnings

from liblibra_core import *
import util.libutil as comn
from libra_py import units
from libra_py import dynamics_plotting
import models
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot

import libra_py.dynamics.exact.compute as dvr
import libra_py.dynamics.exact.save as dvr_save

import libra_py.data_savers as data_savers



#from matplotlib.mlab import griddata
#%matplotlib inline 
warnings.filterwarnings('ignore')

colors = {}
colors.update({"11": "#8b1a0e"})  # red       
colors.update({"12": "#FF4500"})  # orangered 
colors.update({"13": "#B22222"})  # firebrick 
colors.update({"14": "#DC143C"})  # crimson   
colors.update({"21": "#5e9c36"})  # green
colors.update({"22": "#006400"})  # darkgreen  
colors.update({"23": "#228B22"})  # forestgreen
colors.update({"24": "#808000"})  # olive      
colors.update({"31": "#8A2BE2"})  # blueviolet
colors.update({"32": "#00008B"})  # darkblue  
colors.update({"41": "#2F4F4F"})  # darkslategray

clrs_index = ["11", "21", "31", "41", "12", "22", "32", "13","23", "14", "24"]



import argparse
parser = argparse.ArgumentParser(description='Exact...')
parser.add_argument('--model_indx', type=int)
args = parser.parse_args()


##### 1. Choose the model

#################################
# Give the model used an index
model_indx = args.model_indx
################################

model_params = models.all_model_params[model_indx]


##### 2. Choose the Nonadiabatic Dynamics Methodology
NSTATES = model_params["nstates"]

q0 = -2.0
p0 = 10.0
istate = NSTATES-1

prf = F"EXACT-model{model_indx}"

exact_params = { "nsteps":1000, "dt":4.0, "progress_frequency":1.0/10,
                 "rmin":[-30.0], "rmax":[30.0], "dx":[0.025], "nstates":NSTATES,
                  "x0":[q0], "p0":[p0], "istate":[1,  istate], "masses":[1800.0], "k":[0.01],
                  "integrator":"SOFT",
                  "mem_output_level":0, "txt_output_level":0, "txt2_output_level":0, "hdf5_output_level":2, 
                  "properties_to_save":[ "timestep", "time", "Epot_dia", "Ekin_dia", "Etot_dia",
                                         "Epot_adi", "Ekin_adi", "Etot_adi", "norm_dia", "norm_adi",
                                         "pop_dia", "pop_adi", "q_dia", "q_adi", "p_dia", "p_adi" ],
                  "prefix":prf, "prefix2":prf, "use_compression":0, "compression_level":[0, 0, 0]
               }

wfc = dvr.init_wfc(exact_params, models.potential, model_params)
savers = dvr_save.init_tsh_savers(exact_params, model_params, exact_params["nsteps"], wfc)
dvr.run_dynamics(wfc, exact_params, model_params, savers)



