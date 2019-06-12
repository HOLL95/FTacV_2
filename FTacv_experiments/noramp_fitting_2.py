import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_noramp  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]
dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Black"#"/Gold/Large"#
Method ="N_"#"GoldLarge_1"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        results2=np.loadtxt(path+"/"+data)
dec_amount=8
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=results2[0::dec_amount, 1]

current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=results2[0::dec_amount, 1]

length_list=[1e4, 2e4, 3e4]
dec_list=[8, 16, 32, 64]
repeat_num=5
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 382,  #     (uncompensated resistance ohms)
    'Cdl': 5.440677575193328e-05, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1.9e-10,          # (surface coverage per unit area)
    'k_0': 357, #(reaction rate s-1)
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*math.pi/2,
    "time_end": None,
    'num_peaks':50
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "test": False,
    "likelihood":likelihood_options[0],
    "numerical_method": solver_list[1],
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,9,1),
    "experiment_time": time_results,
    "experiment_current": current_results,
    "experiment_voltage":voltage_results,
    "bounds_val":20,
    "signal_length":int(2e4),
}
param_list['E_0']=0.217740555023939#(param_list['E_reverse']-param_list['E_start'])/2
noramp_fit=single_electron(param_list, simulation_options, other_values)
harm_class=harmonics(other_values["harmonic_range"], noramp_fit.nd_param.omega*noramp_fit.nd_param.c_T0, 0.5)
time_results=noramp_fit.other_values["experiment_time"]
current_results=noramp_fit.other_values["experiment_current"]
voltage_results=noramp_fit.other_values["experiment_voltage"]
noramp_fit.label="MCMC"
noramp_fit.optim_list=["k_0", "omega"]
num_params=5

param_ranges={
"E_0_range":[0.2, 0.3, 0.4, 0.5, 0.6],#np.linspace(estart+(de/2), ereverse-(de/2), num_params),
"k_0_range":[0.2,2, 20, 200, 2000],
"Ru_range":[1, 10, 100, 1000, 5000],#[10,100,500,1000,1500,2000,2500,3000],
"Cdl_range":[0,1e-7, 1e-6, 1e-5, 1e-4],
"CdlE1_range":[-2, -0.5,0, 0.5, 2],
"CdlE2_range":[-0.05, -0.01,0, 0.01, 0.05],
"CdlE3_range":[-0.01, -0.005,0, 0.005, 0.01],
"alpha_range":[0.1, 0.2,0.5, 0.7, 0.9],
}
param_keys=param_ranges.keys()
true_params=['E_0', 'k_0', 'Ru','Cdl', 'CdlE1','CdlE2','CdlE3','alpha']
start=time.time()
end=np.where(noramp_fit.time_vec>0.01)
noramp_fit.time_vec=np.array(noramp_fit.time_vec)
cdl_count=0
voltages=noramp_fit.define_voltages()
voltages=voltages
sort_idx=np.argsort(voltages)
voltages=np.multiply(voltages, noramp_fit.nd_param.c_E0)
plt.plot(voltages)
plt.plot(np.multiply(voltage_results, noramp_fit.nd_param.c_E0))
plt.show()
voltages=voltages[end]
num_times=5
means=[0.217740555023939, 357.3609131669447, 0, 5.440677575193328e-05, 0, 1.9319598326712751e-11, 8.94098189688349, 3*math.pi/2, 0.9000000005558953]
means2=[0.217740555023939, 357.3609131669447, 0, 5.440677575193328e-05, 0, 1.9319598326712751e-11, 8.94098189688349, 3*math.pi/2, 0.9000000005558953]
list1=['E_0', 'k_0', 'Ru','Cdl', 'CdlE1','gamma','omega', 'phase','alpha']
noramp_fit.optim_list=list1
time_series=noramp_fit.test_vals(means, "timeseries", test=False)
noramp_fit.simulation_options["numerical_method"]="inverted"
time_series2=noramp_fit.test_vals(means2, "timeseries", test=False)
time_series=time_series
#time_series=np.flip(time_series*-1)
voltage_results=np.multiply(voltage_results, noramp_fit.nd_param.c_E0)
plt.plot(voltage_results[end], np.multiply(time_series[end], noramp_fit.nd_param.c_I0), label="flipped")
plt.plot(voltage_results[end], (np.multiply(time_series2[end], noramp_fit.nd_param.c_I0)), label="original")
#plt.plot(voltage_results[end], np.flip(np.multiply(time_series[end], noramp_fit.nd_param.c_I0)))
plt.legend()
plt.show()
