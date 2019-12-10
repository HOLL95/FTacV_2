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
dec_amount=32
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
sf_list=[200, 400, 800, 1000, 2000, 5000, 10000, 20000, 50000]
prev_data="data"

omega=10
sf=100
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
'num_peaks': 2
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
syn_fit=single_electron(None, param_list, simulation_options, other_values, param_bounds)
true_vals=[0.25, 100, 300, 1e-10, 1e-5, 0,0, 3*math.pi/2, 3*math.pi/2, 0.5]
syn_fit.def_optim_list(["E_0", "k_0", "Ru", "gamma", "Cdl", "CdlE1", "CdlE2", "phase", "cap_phase", "alpha"])
data2=syn_fit.test_vals(true_vals, "timeseries")
voltages2=syn_fit.define_voltages()
times2=syn_fit.time_vec
after_transient=5/param_list["omega"]
time_idx2=tuple(range(0, len(data2)))#tuple(np.where(times2>after_transient))
#plt.subplot(1,2,1)
#plt.plot(syn_fit.e_nondim(voltages[time_idx1]), syn_fit.i_nondim(data[time_idx1]), label=("v "+str(sf)))
plot_voltages=(voltages2)
plot_current=(data2)

plot_times=syn_fit.t_nondim(syn_fit.time_vec)
variables=["Cdl", "dE", "Ru", "d_I", "gamma", "k_0", "theta", "alpha", "E", "I","E_0"]
sensitive_params=["E_0","Ru", "Cdl", "k_0", "alpha"]
symbolic_vars=[sym.symbols(x) for x in variables]
symbol_dict=dict(zip(variables, symbolic_vars))
d_theta_1=((1-symbol_dict["theta"])*symbol_dict["k_0"]*sym.exp((1-symbol_dict["alpha"])*(symbol_dict["E"]-(symbol_dict["Ru"]*symbol_dict["I"])-symbol_dict["E_0"])))
d_theta_2=((symbol_dict["theta"])*symbol_dict["k_0"]*sym.exp(-symbol_dict["alpha"]*(symbol_dict["E"]-(symbol_dict["Ru"]*symbol_dict["I"])-symbol_dict["E_0"])))
d_theta=(d_theta_1-d_theta_2)
cap=symbol_dict["Cdl"]*(symbol_dict["dE"]-(symbol_dict["Ru"]*symbol_dict["d_I"]))
F=cap+(symbol_dict["gamma"]*d_theta)
G=d_theta
lambda_F_vars=[str(x) for x in list(F.free_symbols)]
lambda_F=sym.lambdify(lambda_F_vars, F)
state_variables=[symbol_dict["I"], symbol_dict["theta"]]
J_F=[sym.diff(F, x) for x in state_variables]
J_G=[sym.diff(G, x) for x in state_variables]
dF_dp=[sym.diff(F, x) for x in sensitive_params]
dG_dp=[sym.diff(G, x) for x in sensitive_params]
init_sensitivity=[0, 0]
I_sense=np.zeros((len(sensitive_params), len(voltages2)))
theta_sense=np.zeros((len(sensitive_params), len(voltages2)))
func_lists=[J_F, J_G, dF_dp, dG_dp]
keys=["J_F", "J_G", "dF_dp", "dG_dp"]
lambda_dict={}
variable_dict={}
print(dF_dp, "**"*50)
def DF_DI(function=lambda_F, **variables, ):
    delta=1e-5
    orig_I=copy.deepcopy(variables["I"])
    variables["I"]=orig_I-(delta*orig_I)
    delta_I_1=function(**variables)#
    variables["I"]=orig_I+(delta*orig_I)
    delta_I_2=function(**variables)
    return (delta_I_2-delta_I_1)/(2*delta*orig_I)
for i in range(0, len(func_lists)):
    derivatives=func_lists[i]
    expr_variables=[[str(x) for x in list(expr.free_symbols)] for expr in derivatives]
    if keys[i][0].lower()=="j":
        params=[str(x) for x in state_variables]
    elif keys[i][0].lower()=="d":
        params=sensitive_params
    lambda_dict[keys[i]]=dict(zip(params, [sym.lambdify(variables, derivative) for variables, derivative in zip(expr_variables, derivatives)]))
    variable_dict[keys[i]]=dict(zip(params, expr_variables))
lambda_dict["J_F"]["I"]=DF_DI
variable_dict["J_F"]["I"]=lambda_F_vars
d_I=np.zeros(len(plot_current))
d_E=np.zeros(len(plot_current))
d_E2=np.zeros(len(plot_current))
for i in range(0, len(plot_current)-1):
    d_I[i]=(plot_current[i+1]-plot_current[i])/syn_fit.nd_param.sampling_freq#(syn_fit.time_vec[i+1]-syn_fit.time_vec[i])
    d_E[i]=isolver_martin_brent.dEdt(syn_fit.nd_param.nd_omega, syn_fit.nd_param.phase, syn_fit.nd_param.d_E, syn_fit.time_vec[i])
theta=syn_fit.calc_theta(plot_current)
def var_dict(desired_params, param_dict):
    return dict(zip(desired_params, [param_dict[x] for x in desired_params]))

for i in range(1, len(voltages2)):
    vals=[syn_fit.nd_param.Cdl, d_E[i], syn_fit.nd_param.Ru, d_I[i], syn_fit.nd_param.gamma, syn_fit.nd_param.k_0, theta[i], syn_fit.nd_param.alpha, plot_voltages[i], plot_current[i], syn_fit.nd_param.E_0]
    val_dict=dict(zip(variables, vals))
    d_s=np.zeros(len(sensitive_params))
    jacobian_entry_I=lambda_dict["J_F"]["I"](**var_dict(variable_dict["J_F"]["I"], val_dict))# NEED TO MULTIPLY BY CURRENT SENSITIVITY MORON
    jacobian_entry_theta=lambda_dict["J_F"]["theta"](**var_dict(variable_dict["J_F"]["theta"], val_dict))
    for j in range(0, len(sensitive_params)):
        d_func_entry=dF_dp[j].evalf(subs=val_dict)
        d_s[j]=d_func_entry
    I_sense[:,i]=d_func_entry#np.add(I_sense[:,i-1],np.multiply(syn_fit.nd_param.sampling_freq, d_s))

for i in range(0, len(sensitive_params)):
    plt.subplot(len(sensitive_params), 1, i+1)
    plt.plot(I_sense[i, :])
plt.show()
"""
#didt=((symbol_dict["Cdl"]*symbol_dict["dE"])+(symbol_dict["gamma"]*d_theta))/(symbol_dict["Cdl"]*symbol_dict["Ru"])
sym.pprint(didt)
derivatives=[sym.diff(didt, x) for x in sensitive_params]
#derivatives[0]=derivatives[1]
[sym.pprint(deriv) for param, deriv in zip(sensitive_params, derivatives)]


print(expr_variables)
lambdas=[sym.lambdify(variables, derivative) for variables, derivative in zip(expr_variables, derivatives)]
sensitivity_mat=np.zeros((len(times2), len(derivatives)))
theta=syn_fit.calc_theta(plot_current)

plt.plot(np.subtract(d_E, d_E2))
plt.title("numerical_derivs")
plt.show()
time_idx=tuple(np.where(syn_fit.t_nondim(syn_fit.time_vec)>(2/omega)))
"""
