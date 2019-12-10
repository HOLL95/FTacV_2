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
dec_amount=32
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
sf_list=[200, 400, 800, 1000, 2000, 5000, 10000, 20000, 50000]
prev_data="data"

omega=10
sf=2000
param_list={
"E_0":0.2,
'E_start': estart, #(starting dc voltage - V)
'E_reverse': ereverse,
'omega':omega,#8.88480830076,  #    (frequency Hz)
"original_omega":omega,
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
'num_peaks': 10
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
time_start=5/(param_list["omega"])
#start_idx=np.where(time_results1>time_start)
simulation_options={
"no_transient":1/omega,
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
"bounds_val":25,
"signal_length":int(2e4)
}
param_bounds={
'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
'omega':[0.95*param_list['omega'],1.05*param_list['omega']],#8.88480830076,  #    (frequency Hz)
'Ru': [0, 1e4],  #     (uncompensated resistance ohms)
'Cdl': [0,1e-4], #(capacitance parameters)
'CdlE1': [-0.05,0.15],#0.000653657774506,
'CdlE2': [-0.01,0.01],#0.000245772700637,
'CdlE3': [-0.01,0.01],#1.10053945995e-06,
'gamma': [1e-11,1e-9],
'k_0': [10, 1e3], #(reaction rate s-1)
'alpha': [0.4, 0.6],
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
syn_fit2=single_electron(None, param_list, simulation_options, other_values, param_bounds)
true_vals=[0.25, 200, 300, 1e-10, 1e-5, -1e-4, 1e-4, 3*math.pi/2, 3*math.pi/2, 0.5]
syn_fit2.def_optim_list(["E_0", "k_0", "Ru", "gamma", "Cdl", "CdlE1", "CdlE2", "phase", "cap_phase", "alpha"])
data2=syn_fit2.test_vals(true_vals, "timeseries")
voltages2=syn_fit2.define_voltages()

times2=syn_fit2.time_vec
after_transient=5/param_list["omega"]
time_idx2=tuple(range(0, len(data2)))#tuple(np.where(times2>after_transient))
#plt.subplot(1,2,1)
#plt.plot(syn_fit2.e_nondim(voltages[time_idx1]), syn_fit2.i_nondim(data[time_idx1]), label=("v "+str(sf)))
plot_voltages=syn_fit2.e_nondim(voltages2)[syn_fit2.time_idx:]
plot_current=syn_fit2.i_nondim(data2)
data=plot_current
plot_times=syn_fit2.t_nondim(syn_fit2.time_vec)
variables=["Cdl", "dE", "Ru", "dI", "gamma", "k_0", "theta", "alpha", "E", "I","E_0"]
noise_pc=0.05
error=max(abs(data2))*noise_pc
noise=np.random.normal(0, error, len(data))
noisy_data=np.add(noise, data2)
#error=np.std(np.abs(np.subtract(noisy_data, data)))

mcmc_problem=pints.SingleOutputProblem(syn_fit2, syn_fit2.time_vec[syn_fit2.time_idx:], noisy_data)
syn_fit2.secret_data_time_series=noisy_data
plt.plot(plot_voltages, noisy_data)
plt.show()
#updated_lb=np.append([noramp_results.param_bounds[key][0] for key in master_optim_list],0.75*error)
#updated_ub=np.append([noramp_results.param_bounds[key][1] for key in master_optim_list], 1.25*error)
updated_lb=np.append(np.multiply(true_vals, 0.5),0.01*error)
updated_ub=np.append(np.multiply(true_vals, 1.5), 10*error)
updated_b=[updated_lb, updated_ub]
updated_b=np.sort(updated_b, axis=0)
log_liklihood=pints.GaussianLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(true_vals, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
syn_fit2.simulation_options["label"]="MCMC"
mcmc = pints.MCMCController(log_posterior, 3, xs,method=pints.HaarioBardenetACMC)
mcmc.set_parallel(False)
syn_fit2.simulation_options["test"]=False
mcmc.set_max_iterations(10000)
chains=mcmc.run()
f=open("Synthetic_MCMC_"+"lower_info"+"_"+"run1", "wb")
np.save(f, chains)
f.close()
pints.plot.trace(chains)
#print file
plt.show()
