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

length_list=[1e4, 2e4, 3e4]
dec_list=[8, 16, 32, 64]
repeat_num=5
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 0.0, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,          # (surface coverage per unit area)
    'k_0': 10, #(reaction rate s-1)
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
    'time_end':1000,
    'num_peaks': 2
}
param_list['E_0']=(param_list['E_reverse']-param_list['E_start'])/2
harmonic_range=range(1,12,1)
noramp_params=single_electron(param_list, params_for_opt, harmonic_range, 0.5)
noramp_params.nd_param.time_end=(noramp_params.nd_param.num_peaks/noramp_params.nd_param.nd_omega)*2*math.pi
noramp_params.numerical_method="Brent minimisation"
noramp_params.bounds_val=10000000
signal_length=int(1e4)
noramp_params.times(signal_length)
frequencies=np.fft.fftfreq(signal_length, noramp_params.time_vec[1]-noramp_params.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
noramp_params.frequencies=frequencies
last_point= (harmonic_range[-1]*noramp_params.nd_param.omega)+(noramp_params.nd_param.omega*0.5)
plot_frequencies=frequencies[np.where(frequencies<last_point)]
noramp_params.test_frequencies=plot_frequencies
harm_class=harmonics(harmonic_range, noramp_params.nd_param.omega*noramp_params.nd_param.c_T0, 0.5)
noramp_params.label="MCMC"
noramp_params.optim_list=["k_0", "omega"]
num_params=4

param_ranges={
"E_0_range":np.linspace(estart+(de/2), ereverse-(de/2), num_params),
"k_0_range":[0.2,2, 20, 200],
"Ru_range":[0, 1.0, 50, 1000],
"Cdl_range":[0, 1e-6, 1e-5, 1e-4],
"CdlE1_range":[0, 1e-6, 1e-5, 1e-4],
"CdlE3_range":[0, 1e-6, 1e-5, 1e-4],
"alpha_range":np.linspace(0.3, 0.7, num_params),
"phase_range":np.linspace(0, 2*math.pi, num_params)
}
param_keys=param_ranges.keys()
true_params=param_list.keys()
start=time.time()
fig, ax=plt.subplots(2,4)
for i in range(0, len(param_keys)):
    for j in range(0, len(true_params)):
        if true_params[j] in param_keys[i]:
            parameter_name=true_params[j].strip()
            break
    noramp_params.optim_list=[parameter_name, "omega"]
    for k in range(0, num_params):
        parameter_val=param_ranges[param_keys[i]][k]
        print parameter_val
        time_series=noramp_params.simulate([parameter_val, 8.94], noramp_params.time_vec, "no", "timeseries", "no")
        check_dict=vars(noramp_params.nd_param)
        print check_dict[parameter_name]
        #ax[i].plot(noramp_params.time_vec, time_series)
    check_dict=vars(noramp_params.nd_param)
    print parameter_name, check_dict[parameter_name]
    noramp_params.nd_param.non_dimensionalise(parameter_name, param_list[parameter_name])
    check_dict=vars(noramp_params.nd_param)
    print parameter_name, check_dict[parameter_name]
    #time_series=noramp_params.simulate([k0_range[i], 8.94], noramp_params.time_vec, "no", "timeseries", "no")
    #plt.plot(noramp_params.time_vec[:], time_series)
