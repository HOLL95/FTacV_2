import numpy as np
import matplotlib.pyplot as plt
import math
import time
start=time.time()
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
import matplotlib.gridspec as gridspec
print("import", time.time()-start)
types=["current", "voltage"]
start=time.time()
noramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={}, dec_amount=4)
ramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={}, dec_amount=64)
print("read", time.time()-start)
val=500*0
simulation_options={
    "no_transient":val,
    "numerical_debugging": False,
    "experimental_fitting":False,
    "dispersion":False,
    "dispersion_bins":40,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"sinusoidal",
    "label": "MCMC",
    "optim_list":[]
}
ramped_simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":False,
    "dispersion":False,
    "dispersion_bins":40,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"ramped",
    "label": "MCMC",
    "optim_list":[]
}
noramp_other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(3,7,1)),
    "experiment_time": None,#noramp_startup.time_results["GC4_1_cv"],
    "experiment_current":None,#noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":None,#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":int(2e5),
}

ramped_other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(1,9,1)),
    "experiment_time": None,#ramp_startup.time_results["GC4_1_ramp_cv"],
    "experiment_current": None,#ramp_startup.current_results["GC4_1_ramp_cv"],
    "experiment_voltage":None,#ramp_startup.voltage_results["GC4_1_ramp_cv"],
    "bounds_val":20,
    "signal_length":int(5e5),
}
param_bounds={
    'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.95*noramp_startup.generic_noramp_params['omega'],1.05*noramp_startup.generic_noramp_params['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [40, 700],  #     (uncompensated resistance ohms)
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
#plt.rcParams.update({'font.size': 12})
noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
#noramp_simulations=single_electron(ramp_startup.generic_ramped_params, ramped_simulation_options, ramped_other_values)
harm_class=harmonics(noramp_other_values["harmonic_range"],noramp_simulations.nd_param.omega, 0.2)
noramp_simulations.def_optim_list(["k_0", "omega"])
omegas=[9, 9*25]
k_s=[5, 10, 50, 100, 150, 300, 600]
def RMSE(series1, series2):
    return np.sqrt((np.sum(1/(len(series1))*np.power(np.subtract(series1, series2),2))))
def relative_error(current, previous):
    sub=np.subtract(current, previous)
    div=np.divide(sub, current)
    return np.sum(abs(div))
s_max=[200, 2000]
s_interval=[1, 10]
fig, ax=plt.subplots(1,2)
for q in range(0, 2):
    sampling_freq=range(20,s_max[q] ,s_interval[q])
    sfs=[1/x for x in sampling_freq]
    print(sfs)
    #
    #k_s=np.linspace(0.1, 0.3, 5)
    fig = plt.figure(figsize=(9,9))


    colors=["green", "gold", "red", "gold"]
    alignment=["bottom", "left", "right", "top"]
    outer_val=0.005
    inner_val=0.05
    plot_axes=[]
    noramp_startup.generic_noramp_params["original_omega"]=8.94
    error_array=np.zeros(len(sfs))
    for i in range(0, len(k_s)):
        noramp_startup.generic_noramp_params["sampling_freq"]=sfs[0]
        noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
        noramp_simulations.def_optim_list(["k_0","omega"])
        previous_series=noramp_simulations.i_nondim(noramp_simulations.test_vals([k_s[i],8.94], "timeseries"))
        previous_time=noramp_simulations.time_vec
        for j in range(1, len(sfs)):
            noramp_startup.generic_noramp_params["sampling_freq"]=sfs[j]
            noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
            noramp_simulations.def_optim_list(["k_0","omega"])

            current_series=noramp_simulations.i_nondim(noramp_simulations.test_vals([k_s[i],8.94], "timeseries"))
            current_time=noramp_simulations.time_vec
            interped_time_series=np.interp(previous_time, current_time, current_series)


            #plt.plot(previous_time, previous_series)
            #plt.plot(previous_time, interped_time_series)
            #plt.show()
            error_array[j]=relative_error(interped_time_series[1:], previous_series[1:])
            previous_series=current_series
            previous_time=current_time
        ax[q].plot(sampling_freq[5:], error_array[5:], label=k_s[i])
        ax[q].set_xlabel("Sampling frequency")
        ax[q].set_ylabel("Relative error")
        ax[q].legend()
plt.show()
