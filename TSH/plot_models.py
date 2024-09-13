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
from libra_py import dynamics_plotting

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



# New plotting:
# Common setups
plot_params = {"figsize":[20, 5], "titlesize":28, "labelsize":20, "fontsize": 18, "xticksize":26, "yticksize":26,
               "colors": colors, "clrs_index": clrs_index,
               "prefix":F"case", "save_figures":1, "do_show":0,
               "xlim":[-4, 4], "ylim":[-0.08, 0.08], "ylim2":[-5, 5], "show_nac_abs":1,
               "plotting_option":1, "nac_idof":0, "margin_bottom":0.2, "margin_left":0.1,
               "linewidth":8 }

# Esch-Levine
#for indx, prms in enumerate([models.all_model_params[6]]):
for indx, prms in enumerate(models.all_model_params):
    list_states = [x for x in range(prms["nstates"])]
    plot_params.update( { "prefix":F"HAM_{indx}" })
    if indx in [5,6]:
        plot_params.update( { "xlim":[-6, 10], "ylim":[-0.005, 0.01], "ylim2":[-5, 5] } ) 
    dynamics_plotting.plot_surfaces(models.compute_model, [ prms ], list_states, -15.0, 10.0, 0.025, plot_params)

