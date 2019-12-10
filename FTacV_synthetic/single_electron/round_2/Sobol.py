import numpy as np
import matplotlib.pyplot as plt
import math
import time
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
dirs=(dir_path.split("/"))
desired_folder="FTacv_experiments"
global_folder="FTacV_2"
desired_path=("/").join(np.append(dirs[1:dirs.index(global_folder)+1], desired_folder))
sys.path.append("/"+desired_path+"/")
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
import isolver_martin_brent
from harmonics_plotter import harmonics
import os
import pickle
import pints
import pints.plot
import sympy as sym
import copy
from SALib.analyze import sobol
from SALib.sample import saltelli
dec_amount=32
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
sf_list=[200, 400, 800, 1000, 2000, 5000, 10000, 20000, 50000]
prev_data="data"

omega=10
sf=1000
param_list={
"E_0":0.2,
'E_start': estart, #(starting dc voltage - V)
'E_reverse': ereverse,
'omega':omega,#8.88480830076,  #    (frequency Hz)
'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
'v': 1e-3,   #       (scan rate s^-1)
'area': 0.07, #(electrode surface area cm^2)
'Ru': 1.0,  #     (uncompensated resistance ohms)
'Cdl': 1e-6, #(capacitance parameters)
'CdlE1': 0,#0.000653657774506,
'CdlE2': 0,#0.000245772700637,
'CdlE3': 0,#1.10053945995e-06,
'gamma': 1e-10,
"original_gamma":1e-10,        # (surface coverage per unit area)
'k_0': 10, #(reaction rate s-1)
'alpha': 0.5,
"E0_mean":0.2,
"E0_std": 0.09,
"k0_shape":0.954,
"k0_loc":100,
"k0_scale":50,
"k0_range":1e3,
"cap_phase":0,
'sampling_freq' : (1.0/sf),
'phase' : 3*(math.pi/2),
"time_end": None,
'num_peaks': 8
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
time_start=5/(param_list["omega"])
#start_idx=np.where(time_results1>time_start)
simulation_options={
"no_transient":False,
"numerical_debugging": False,
"experimental_fitting":False,
"dispersion":False,
"dispersion_bins":16,
"test": False,
"method": "sinusoidal",
"likelihood":likelihood_options[0],
"numerical_method": solver_list[1],
"label": "MCMC",
"optim_list":[]
}
other_values={
"filter_val": 0.5,
"harmonic_range":list(range(3,9,1)),
"experiment_time": None,
"experiment_current": None,
"experiment_voltage":None,
"bounds_val":20,
"signal_length":int(2e4)
}
param_bounds={
'E_0':[estart, ereverse],#[param_list['E_start'],param_list['E_reverse']],
'omega':[0.95*param_list['omega'],1.05*param_list['omega']],#8.88480830076,  #    (frequency Hz)
'Ru': [0, 1e3],  #     (uncompensated resistance ohms)
'Cdl': [0,1e-4], #(capacitance parameters)
'CdlE1': [-0.15,0.15],#0.000653657774506,
'CdlE2': [-0.01,0.01],#0.000245772700637,
'CdlE3': [-0.01,0.01],#1.10053945995e-06,
'gamma': [1e-11,1e-9],
'k_0': [0, 1e4], #(reaction rate s-1)
'alpha': [0.2, 0.8],
"cap_phase":[0, 2*math.pi],
"E0_mean":[0.19, 0.208],
"E0_std": [0.001, 0.2],
"k0_shape":[0,2],
"k0_loc":[0, 1e3],
"k0_scale":[0,1e3],
"k0_range":[1e2, 1e4],
'phase' : [0, 2*math.pi]
}
vs=[1.0/omega, 2*math.pi/omega, omega]

param_list["v"]=1.0/(omega)
param_list["sampling_freq"]=1.0/sf
syn_fit=single_electron(None, param_list, simulation_options, other_values, param_bounds)
true_vals=[0.25, 100, 300, 1e-10, 1e-5, 0,0, 3*math.pi/2, 3*math.pi/2, 0.5]
orig_list=["E_0", "k_0", "Ru", "gamma", "Cdl","CdlE1", "CdlE2", "phase", "cap_phase", "alpha"]
syn_fit.def_optim_list(["E_0", "k_0", "Ru", "gamma", "Cdl","CdlE1", "CdlE2", "phase", "cap_phase", "alpha"])
data2=syn_fit.test_vals(true_vals, "timeseries")
voltages2=syn_fit.define_voltages()
times2=syn_fit.time_vec
after_transient=5/param_list["omega"]
desired_bounds=["E_0", "k_0", "Ru", "Cdl", "alpha"]

for i in range(0, len(desired_bounds)):
    idx=orig_list.index(desired_bounds[i])
    param_bounds[desired_bounds[i]]=[0, 2*true_vals[idx]]

syn_fit.param_bounds=param_bounds
syn_fit.def_optim_list(desired_bounds)
print(syn_fit.param_bounds)
#plt.subplot(1,2,1)
#plt.plot(syn_fit.e_nondim(voltages[time_idx1]), syn_fit.i_nondim(data[time_idx1]), label=("v "+str(sf)))
plot_voltages=(voltages2)
plot_current=(data2)
syn_fit.simulation_options["likelihood"]="timeseries"
syn_fit.simulation_options["label"]="cmaes"
problem={
    "num_vars":len(syn_fit.optim_list),
    "names":syn_fit.optim_list,
    "bounds":[[0, 1]]*len(syn_fit.optim_list)
}
param_vals = saltelli.sample(problem, 10000)
print(len(param_vals))
Y = np.zeros([param_vals.shape[0]])
def RMSE(vals, data):
    sub=np.subtract(vals, data)
    frac=np.sum(np.power(sub, 2))/len(vals)
    return np.sqrt(frac)
for i in range(0, len(param_vals)):
    #print(list(param_vals[i]))
    time_series=syn_fit.simulate(param_vals[i], [])
    Y[i]=RMSE(time_series, data2)
print("pc_change")
Si=sobol.analyze(problem, Y, print_to_console=True)
"""
Parameter S1 S1_conf ST ST_conf
E_0 0.806649 0.021917 0.966729 0.023029
k_0 0.007711 0.008276 0.134534 0.007202
Ru 0.008671 0.004887 0.033040 0.001965
Cdl 0.013382 0.004664 0.023982 0.000987
alpha 0.005703 0.007530 0.056405 0.002582

Parameter_1 Parameter_2 S2 S2_conf
E_0 k_0 0.089748 0.031275
E_0 Ru 0.003719 0.031097
E_0 Cdl 0.007014 0.029475
E_0 alpha 0.012031 0.032523
k_0 Ru -0.002814 0.012188
k_0 Cdl -0.001953 0.012281
k_0 alpha -0.005976 0.012915
Ru Cdl 0.002340 0.006837
Ru alpha 0.002696 0.007013
Cdl alpha -0.000381 0.006175
"""
"""
Parameter S1 S1_conf ST ST_conf
E_0 0.041565 0.008126 0.085295 0.002863
k_0 0.000183 0.000992 0.001532 0.000354
Ru 0.018421 0.004643 0.027320 0.001298
Cdl 0.893879 0.022060 0.936642 0.021598
alpha 0.000068 0.000263 0.000085 0.000029

Parameter_1 Parameter_2 S2 S2_conf
E_0 k_0 -0.005520 0.012443
E_0 Ru -0.003171 0.012847
E_0 Cdl 0.036073 0.013075
E_0 alpha -0.005562 0.012450
k_0 Ru 0.000039 0.001483
k_0 Cdl 0.000232 0.001695
k_0 alpha 0.000100 0.001390
Ru Cdl 0.001284 0.008787
Ru alpha 0.000352 0.006551
Cdl alpha -0.001240 0.024208

"""
