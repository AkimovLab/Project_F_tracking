#!/usr/bin/env python
# coding: utf-8


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
import models
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.dynamics.tsh.plot as tsh_dynamics_plot
import libra_py.data_savers as data_savers

from recipes import fssh, fssh2, fssh3, gfsh


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
parser = argparse.ArgumentParser(description='TSH...')
parser.add_argument('--model_indx', type=int)
parser.add_argument('--dt', type=float)
parser.add_argument('--method_indx', type=int)
args = parser.parse_args()



##### 1. Choose the model

#################################
# Give the model used an index
model_indx = args.model_indx
################################

model_params = models.all_model_params[model_indx]


##### 2. Choose the Nonadiabatic Dynamics Methodology 

#################################
# Choose the dt and the number of time-steps
# but so that the total trajectory length is the same
dt = args.dt
NSTEPS = int(1000.0/dt)

if model_indx in [5,6]:
    NSTEPS = int(10000.0/dt)

################################

NSTATES = model_params["nstates"]
dyn_general = { "nsteps":NSTEPS, "ntraj":200, "nstates":NSTATES,
                "dt":dt, "num_electronic_substeps":1, "isNBRA":0, "is_nbra":0,
                "progress_frequency":0.1, "which_adi_states":range(NSTATES), "which_dia_states":range(NSTATES),      
                "mem_output_level":3,
                "properties_to_save":[ "timestep", "time", "q", "p", "f", "Cadi", "Cdia", "Epot_ave", "Ekin_ave", "Etot_ave",
                "states", "se_pop_adi", "se_pop_dia", "sh_pop_adi", "sh_pop_dia"],
                "prefix":"adiabatic_md", "prefix2":"adiabatic_md"
              }

#################################
# Give the recipe above an index
method_indx = args.method_indx
#################################

if method_indx == 0:
    gfsh.load(dyn_general)  # GFSH, LD - this is our references
    dyn_general.update({"state_tracking_algo":-1, "do_phase_correction":1 })  # default
elif method_indx == 1:
    gfsh.load(dyn_general) # GFSH, state tracking + phase corrections
    dyn_general.update({"state_tracking_algo":2, "do_phase_correction":1 })  # standard state tracking
    dyn_general.update({"rep_tdse":1, "electronic_integrator":3 })   # 1-point, Hvib integration, with exp_
elif method_indx == 2:
    gfsh.load(dyn_general) # GFSH, with F-tracking
    dyn_general.update({"state_tracking_algo":4, "do_phase_correction":0 })  # standard state tracking
    dyn_general.update({"rep_tdse":1, "electronic_integrator":3, "do_nac_phase_correction":1 })   # 1-point, Hvib integration, with exp_
elif method_indx == 3:
    gfsh.load(dyn_general) # GFSH, no tracking, no phase correction, using Hvib integration
    dyn_general.update({"state_tracking_algo":0, "do_phase_correction":0 })  # standard state tracking
    dyn_general.update({"rep_tdse":1, "electronic_integrator":3, "do_nac_phase_correction":0 })   # 1-point, Hvib integration, with exp_
elif method_indx == 4:
    gfsh.load(dyn_general) # GFSH, with F-tracking, but no phase correction
    dyn_general.update({"state_tracking_algo":4, "do_phase_correction":0 })  # standard state tracking
    dyn_general.update({"rep_tdse":1, "electronic_integrator":3, "do_nac_phase_correction":0 })   # 1-point, Hvib integration, with exp_



### 3. Choose the initial conditions: Nuclear and Electronic

# Nuclear initial conditions 
nucl_params = { "ndof":1,  "q":[-2.0], "p":[10.0], "mass":[1800.0], "force_constant":[ 0.01], "init_type":3 }

if model_indx==5:
    nucl_params.update( {"q":[7.5], "p":[0.0] } )
elif model_indx==6:
    nucl_params.update( {"q":[-4.0], "p":[0.0] } )

# Electronic initial conditions - we start on the top state 
istate = NSTATES - 1
if model_indx==5:
    istate = 1
elif model_indx==6:
    istate = 0


istates = [0.0 for i in range(NSTATES)]
istates[istate] = 1.0

    
elec_params = {"verbosity":2, "init_dm_type":0,"ndia":NSTATES, "nadi":NSTATES, 
               "rep":1, "init_type":0, "istates":istates, "istate":istate
              }

### 4. Run the calculations and save the results into the HDF5 files

dyn_params = dict(dyn_general)
pref = F"TSH-method{method_indx}-model{model_indx}-dt{dt}"
dyn_params.update({ "prefix":pref, "prefix2":pref })

print(F"Computing {pref}")    

#sys.exit(0)

rnd = Random()
res = tsh_dynamics.generic_recipe(dyn_params, models.compute_model, model_params, elec_params, nucl_params, rnd)


### 5. Plot the results - preliminary plotting
pref = F"TSH-method{method_indx}-model{model_indx}-dt{dt}"

NTRAJ = dyn_general["ntraj"]
plot_params = { "prefix":pref, "filename":"mem_data.hdf", "output_level":3,
                "which_trajectories":list(range(NTRAJ)), "which_dofs":[0], "which_adi_states":list(range(NSTATES)), 
                "which_dia_states":list(range(NSTATES)), 
                "frameon":True, "linewidth":3, "dpi":300,
                "axes_label_fontsize":(8,8), "legend_fontsize":8, "axes_fontsize":(8,8), "title_fontsize":8,
                "what_to_plot":["coordinates", "momenta",  "forces", "energies", "phase_space", "se_pop_adi",
                                "se_pop_dia", "sh_pop_adi", "sh_pop_dia" ], 
                "which_energies":["potential", "kinetic", "total"],
                "save_figures":1, "do_show":0, "no_label":1
              }

#tsh_dynamics_plot.plot_dynamics(plot_params)

