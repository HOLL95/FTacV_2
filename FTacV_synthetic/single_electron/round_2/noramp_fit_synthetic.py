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

omega=150
sf=2000
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
'num_peaks': 10
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
plot_voltages=syn_fit2.e_nondim(voltages2)
plot_current=syn_fit2.i_nondim(data2)
plt.plot(plot_voltages, plot_current)
plt.show()
plot_times=syn_fit2.t_nondim(syn_fit2.time_vec)
variables=["Cdl", "dE", "Ru", "dI", "gamma", "k_0", "theta", "alpha", "E", "I","E_0"]
sensitive_params=["Cdl", "Ru", "k_0", "alpha", "E_0"]
symbolic_vars=[sym.Symbol(x) for x in variables]
symbol_dict=dict(zip(variables, symbolic_vars))
d_theta_1=(1-symbol_dict["theta"])*symbol_dict["k_0"]*sym.exp((1-symbol_dict["alpha"])*(symbol_dict["E"]-symbol_dict["Ru"]*symbol_dict["I"]-symbol_dict["E_0"]))
d_theta_2=(symbol_dict["theta"])*symbol_dict["k_0"]*sym.exp(-symbol_dict["alpha"]*(symbol_dict["E"]-symbol_dict["Ru"]*symbol_dict["I"]-symbol_dict["E_0"]))
d_theta=d_theta_1+d_theta_2
cap=symbol_dict["Cdl"]*(symbol_dict["dE"]-symbol_dict["Ru"]*symbol_dict["dI"])
didt=cap+(symbol_dict["gamma"]*d_theta)
print(didt)
derivatives=[sym.diff(didt, x) for x in sensitive_params]
expr_variables=[list(expr.free_symbols) for expr in derivatives]
lambdas=[sym.lambdify(variables, derivative) for variables, derivative in zip(expr_variables, derivatives)]
print(lambdas)
sensitivity_mat=np.zeros((len(times2), len(derivatives)))
theta=syn_fit2.calc_theta(plot_current)
d_I=np.zeros(len(plot_current))
d_E=np.zeros(len(plot_current))
for i in range(1, len(plot_current)):
    d_I[i-1]=(plot_current[i]-plot_current[i-1])/(syn_fit2.time_vec[i]-syn_fit2.time_vec[i-1])
    d_E[i-1]=isolver_martin_brent.dEdt(syn_fit2.nd_param.nd_omega, syn_fit2.nd_param.phase, syn_fit2.nd_param.d_E, syn_fit2.time_vec[i])
print(len(times2))


for i in range(0, len(times2)-1):
    vals=[syn_fit2.nd_param.Cdl, d_E[i], syn_fit2.nd_param.Ru, d_I[i], syn_fit2.nd_param.gamma, syn_fit2.nd_param.k_0, theta[i], syn_fit2.nd_param.alpha, voltages2[i], plot_current[i], syn_fit2.nd_param.E_0]
    val_dict=dict(zip(variables, vals))
    start=time.time()
    sensitivities=[function(*[val_dict[str(x)] for x in var_i]) for var_i, function in zip(expr_variables, lambdas)]
    sensitivity_mat[i, :]=sensitivities
time_idx=tuple(np.where(syn_fit2.time_vec>(2/omega)))
for i in range(0, len(sensitive_params)):
    plt.subplot(2,3, i+1)
    ax=plt.gca()
    ax.plot(voltages2[time_idx], plot_current[time_idx])
    ax2=ax.twinx()
    ax.set_title(sensitive_params[i])
    ax2.plot(voltages2[time_idx], sensitivity_mat[:, i][time_idx], color="red", linestyle="--")
plt.show()


noise_pc=0.01
error=max(abs(data))*noise_pc
noise=np.random.normal(0, max(abs(data))*noise_pc, len(data))
noisy_data=np.add(noise, data)
#error=np.std(np.abs(np.subtract(noisy_data, data)))
mcmc_problem=pints.SingleOutputProblem(syn_fit2, syn_fit2.time_vec, noisy_data)

#updated_lb=np.append([noramp_results.param_bounds[key][0] for key in master_optim_list],0.75*error)
#updated_ub=np.append([noramp_results.param_bounds[key][1] for key in master_optim_list], 1.25*error)
updated_lb=np.append(np.multiply(true_vals, 0.5),0.5*error)
updated_ub=np.append(np.multiply(true_vals, 1.5), 1.5*error)
updated_b=[updated_lb, updated_ub]
updated_b=np.sort(updated_b, axis=0)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(true_vals, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
syn_fit2.simulation_options["label"]="MCMC"
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_parallel(True)
mcmc.set_max_iterations(20000)
chains=mcmc.run()
f=open("Synthetic_MCMC_"+"lower_info"+"_"+"run1", "wb")
np.save(f, chains)
f.close()
pints.plot.trace(chains)
#print file
plt.show()
