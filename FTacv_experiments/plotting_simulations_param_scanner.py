import numpy as np
import matplotlib.pyplot as plt
import math
import time
start=time.time()
from single_e_class_unified import single_electron, paralell_class
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
import isolver_martin_brent
print("import", time.time()-start)
types=["current", "voltage"]
start=time.time()
noramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={"GC4_1_cv":types, "GC4_2_cv":types, "GC4_3_cv":types}, dec_amount=4)
ramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={"GC4_1_ramp_cv":types,}, dec_amount=64)
print("read", time.time()-start)
simulation_options={
    "no_transient":1/5,
    "numerical_debugging": False,
    "experimental_fitting":False,
    "dispersion":False,
    "dispersion_bins":20,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"sinusoidal",
    "label": "MCMC",
    "optim_list":[]
}

noramp_other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(3,7,1)),
    "experiment_time": False,#noramp_startup.time_results["GC4_1_cv"],
    "experiment_current":False, #noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":False,#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":int(5e3),
}

noramp_startup.generic_noramp_params['omega']=5
noramp_startup.generic_noramp_params['Cdl']=1e-5
noramp_startup.generic_noramp_params['num_peaks']=20
noramp_startup.generic_noramp_params['sampling_freq']=(1.0/400)
param_bounds={
    'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.1*noramp_startup.generic_noramp_params['omega'],1.9*noramp_startup.generic_noramp_params['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 700],  #     (uncompensated resistance ohms)
    'Cdl': [0,0], #(capacitance parameters)
    'CdlE1': [-0.05,0.15],#0.000653657774506,
    'CdlE2': [-0.01,0.01],#0.000245772700637,
    'CdlE3': [-0.01,0.01],#1.10053945995e-06,
    'gamma': [1e-11,1e-9],
    'k_0': [10, 5e3], #(reaction rate s-1)
    'alpha': [0.35, 0.65],
    "cap_phase":[0, 2*math.pi],
    "E0_mean":[0.15, 0.3],
    "E0_std": [0.001, 0.2],
    "k0_shape":[0,2],
    "k0_loc":[0, 5e3],
    "k0_scale":[0,5e3],
    "k0_range":[1e2, 1e4],
    'phase' : [0, 2*math.pi]
}
unit_dict={
    "E_0": "V",
    'E_start': "V", #(starting dc voltage - V)1254684939476, 1.4485673957633017e-10, 8.94084836566341, 4.898877306283271, 4.379038086713167, 0.6999999953311751]
    'E_reverse': "V",
    'omega':"Hz",#8.88480830076,  #    (frequency Hz)
    'd_E': "V",   #(ac voltage amplitude - V) freq_range[j],#
    'v': r'$s^{-1}$',   #       (scan rate s^-1)
    'area': r'$cm^{2}$', #(electrode surface area cm^2)
    'Ru': r'$\Omega$',  #     (uncompensated resistance ohms)
    'Cdl': "F", #(capacitance parameters)
    'CdlE1': "",#0.000653657774506,
    'CdlE2': "",#0.000245772700637,
    'CdlE3': "",#1.10053945995e-06,
    'gamma': r'$mol$ $cm^{-2}$',
    'k_0': r'$s^{-1}$', #(reaction rate s-1)
    'alpha': "",
    "E0_mean":"V",
    "E0_std": "V",
    "k0_shape":"",
    "k0_loc":"",
    "k0_scale":"",
    "cap_phase":"rads",
    'phase' : "rads",
}
noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
noramp_simulations.def_optim_list(["E_0","k_0","Ru","alpha", "Cdl"])
#test=paralell_class(noramp_simulations.nd_param_dict, noramp_simulations.time_vec, "sinusoidal", noramp_simulations.bounds_val, isolver_martin_brent.brent_current_solver)
#noramp_simulations.current_params=noramp_simulations.nd_param_dict
#test.paralell_dispersion(range(20))
#start=time.time()
#test=noramp_simulations.test_vals([0.2, 100, 10, 0.5, 0], "timeseries")
#print(time.time()-start)
#print(len(test))
#plt.plot(test)
#plt.show()
noramp_simulations.param_scanner(noramp_simulations.optim_list, unit_dict,"Scans", 6)
#harm_class=harmonics(noramp_other_values["harmonic_range"],noramp_simulations.nd_param.omega, 0.2)
