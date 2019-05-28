import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_noramp  import single_electron
from params_class import params
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time

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
    'Cdl': 1e-5, #(capacitance parameters)
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
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson"]
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":False,
    "numerical_method": solver_list[1],
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,12,1),
    "experiment_time": None,
    "experiment_current": None,
    "experiment_voltage":None,
    "bounds_val":10000000,
    "signal_length":int(5e4),
    "dec_val":8
}
param_list['E_0']=(param_list['E_reverse']-param_list['E_start'])/2
noramp_params=single_electron(param_list, simulation_options, other_values)
harm_class=harmonics(other_values["harmonic_range"], noramp_params.nd_param.omega*noramp_params.nd_param.c_T0, 0.5)
noramp_params.optim_list=["k_0", "omega"]
num_params=4
param_ranges={
"E_0_range":[0.1, 0.2, 0.3, 0.4],#np.linspace(estart+(de/2), ereverse-(de/2), num_params),
"k_0_range":[0.2,2, 20, 200],
"Ru_range":[ 1.0, 10,100, 1000],
"Cdl_range":[0, 1e-6, 1e-5, 1e-4],
"CdlE1_range":[0, 1e-3, 1e-2, 1e-1],
"CdlE3_range":[0, 1e-5, 1e-4, 1e-3],
"alpha_range":[0.3, 0.4, 0.5, 0.6],
"phase_range":[0, math.pi/2, math.pi, 3*math.pi/2]
}
param_keys=param_ranges.keys()
true_params=param_list.keys()
start=time.time()
end=np.where(noramp_params.time_vec>0.01)
noramp_params.time_vec=np.array(noramp_params.time_vec)
cdl_count=0
for i in range(0, len(param_keys)):
    for j in range(0, len(true_params)):
        if true_params[j] in param_keys[i]:
            if (true_params[j]=="Cdl") and (cdl_count>0):
                continue
            elif (true_params[j]=="Cdl"):
                cdl_count+=1
            parameter_name=true_params[j].strip()
            break
    noramp_params.optim_list=[parameter_name, "omega"]
    for k in range(0, num_params):
        parameter_val=param_ranges[param_keys[i]][k]
        time_series=noramp_params.simulate([parameter_val, 8.94], noramp_params.time_vec, "no", "timeseries", "no")
        time_series=np.array(time_series)
        check_dict=vars(noramp_params.nd_param)
        plt.subplot(2,4,i+1)
        plt.title(parameter_name)
        if parameter_name != "phase":
            plt.plot(noramp_params.time_vec[end], time_series[end], label=str(parameter_val))
        else:
            label=parameter_val/math.pi
            plt.plot(noramp_params.time_vec[end], time_series[end], label=str(label)+r'$\pi$')
        print param_list["CdlE1"], noramp_params.dim_dict["CdlE1"],noramp_params.nd_param.CdlE1
    check_dict=vars(noramp_params.nd_param)
    plt.legend()
    #noramp_params.nd_param.non_dimensionalise(parameter_name, param_list[parameter_name])

    noramp_params=single_electron(param_list, simulation_options, other_values)

    check_dict=vars(noramp_params.nd_param)
plt.show()
