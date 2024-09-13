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


def init_wfc_custom(_params, _potential, _model_params):
    """
    This is a generic and a bit excessive initialization procedure - some of the parameters
    may not be needed for some of the methods, but it allows us to have a simple structure of the program

    istate = [rep, index_of_the_state]

    """

    critical_params = []
    default_params = { "nsteps":200, "dt":10.0, "progress_frequency":0.1,
                       "rmin":[-15.0], "rmax":[15.0], "dx":[0.1], "nstates":2,
                       "x0":[0.0], "p0":[0.0], "istate":[1,0], "masses":[2000.0], "k":[0.001]
                      }
    comn.check_input(_params, default_params, critical_params)

    # Grid properties
    dx = Py2Cpp_double(_params["dx"])
    rmin = Py2Cpp_double(_params["rmin"])
    rmax = Py2Cpp_double(_params["rmax"])
    nstates = _params["nstates"]

    # Dynamical properties
    nsteps = _params["nsteps"]
    dt = _params["dt"]

    # Properties of the initial wavefunction
    istate = _params["istate"]
    k = _params["k"]
    masses = Py2Cpp_double(_params["masses"])

    x0 = Py2Cpp_double(_params["x0"])
    p0 = Py2Cpp_double(_params["p0"])
    ndof = len(x0)

    dx0_tmp = []
    for idof in range(ndof):
        sigmax2 = 0.5*math.sqrt(1.0/(k[idof]*masses[idof]))
        dx0_tmp.append(math.sqrt(sigmax2))
    dx0 = Py2Cpp_double(dx0_tmp)



    # Here we initialize the grid and wavefunction
    wfc = Wfcgrid2( rmin, rmax,  dx, nstates)

    wfc.direct_allocate_tmp_vars(0)    # last time-step wavefunction, diabatic rep
    wfc.direct_allocate_tmp_vars(1)    # last time-step wavefunction, adiabatic rep


    wfc.update_Hamiltonian(_potential, _model_params, 0)  # update Hamiltonian: diabatic -
                                                        # need to compute diabatic-to-adiabatic transform
    wfc.update_Hamiltonian(_potential, _model_params, 1)  # update Hamiltonian: adiabatic, NACs

    wfc.update_propagator_H(0.5*dt)                     # copute the dia-to-adi transform + exp(-i* V *dt/2)

    sq2 = math.sqrt(0.5)
    wfc.add_wfc_Gau(x0, p0, dx0, 1, (1.0+0.0j)*sq2, 3)   # Add to the diabatic or adiabatic state
    wfc.add_wfc_Gau(x0, p0, dx0, 1,-(1.0+0.0j)*sq2, 1)   # Add to the diabatic or adiabatic state

    #if istate[0]==0:
    #    # If we added a diabatic state:
    #    wfc.update_adiabatic()        # update adiabatic wavefunction

    #elif istate[0]==1:
    # If we added an adiabatic state:
    wfc.update_diabatic()        # update diabatic wavefunction


    wfc.update_propagator_K(dt, masses) # update reci space propagator in diabatic rep, exp(-iTdt)
    wfc.update_reciprocal(1)     # update reci of adiabatic function
    wfc.update_reciprocal(0)     # update reci of diabatic function

    print( "Norm (dia) = ", wfc.norm(0) )
    print( "Norm (adi) = ", wfc.norm(1) )
    print( "Ekin (dia) = ", wfc.e_kin(masses, 0) )
    print( "Ekin (adi) = ", wfc.e_kin(masses, 1) )
    print( "Epot (dia) = ", wfc.e_pot(0) )
    print( "Epot (adi) = ", wfc.e_pot(1) )


    return wfc




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

nsteps = 1000
if model_indx==5:
    q0 = 7.5
    p0 = 0.0
    istate = 1
    nsteps = 10000

elif model_indx==6:
    q0 = -4.0
    p0 = 0.0;
    istate = 0
    nsteps = 10000

prf = F"EXACT-model{model_indx}"

exact_params = { "nsteps":nsteps, "dt":4.0, "progress_frequency":1.0/10,
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



