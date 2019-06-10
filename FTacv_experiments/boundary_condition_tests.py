import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_noramp  import single_electron
from model_functions import electrochem_funcs
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,    #  (reverse dc voltage - V)
    'E_0':0.3,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.1, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-6, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,          # (surface coverage per unit area)
    'k_0': 10, #(reaction rate s-1)
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
    "time_end": None,
    'num_peaks': 2
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":False,
    "test":False,
    "likelihood":likelihood_options[0],
    "numerical_method": solver_list[3],
    "label": "MCMC",
    "optim_list":["Ru"]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(3,9,1),
    "bounds_val":2000,
    "signal_length":int(2e4),
}
noramp_bounds=single_electron(param_list, simulation_options, other_values)
bound_test=electrochem_funcs(noramp_bounds.nd_param)
test=noramp_bounds.simulate([1.0], noramp_bounds.time_vec, normalise=False, likelihood="timeseries", test=False )
voltages=noramp_bounds.define_voltages()
kox_vec=np.zeros(other_values["signal_length"])
kred_vec=np.zeros(other_values["signal_length"])
ru_range=[1, 10, 100, 500, 1000, 2000, 3000, 4000]
for j in range(0, 8):
    plt.subplot(2,4, j+1)
    test=noramp_bounds.simulate([ru_range[j]], noramp_bounds.time_vec, normalise=False, likelihood="timeseries", test=False )
    for i in range(0, other_values["signal_length"]):
        kox_vec[i]=bound_test.kox(voltages[i], test[i])
        kred_vec[i]=bound_test.kred(voltages[i], test[i])
    plt.plot(kox_vec)
    plt.plot(kred_vec)
plt.show()
