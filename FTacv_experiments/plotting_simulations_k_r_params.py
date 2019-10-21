import numpy as np
import matplotlib.pyplot as plt
import math
import time
start=time.time()
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
print("import", time.time()-start)
types=["current", "voltage"]
start=time.time()
noramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={}, dec_amount=4)
ramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={}, dec_amount=64)
print("read", time.time()-start)
val=500
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
    "signal_length":int(1e5),
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
plt.rcParams.update({'font.size': 14})
noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
#noramp_simulations=single_electron(ramp_startup.generic_ramped_params, ramped_simulation_options, ramped_other_values)
harm_class=harmonics(noramp_other_values["harmonic_range"],noramp_simulations.nd_param.omega, 0.2)
noramp_simulations.def_optim_list(["k_0", "omega"])
Rus=[50, 100, 150, 200]
k_s=[75, 100, 125, 150]
#k_s2=[1000, 2000, 4000, 8000, 16000]

random_linestyles=["-.", "-", "--", ":"]
random_colours=["b","g", "r", "c" ]
#k_s=np.linspace(0.1, 0.3, 5)
#fig, ax =plt.subplots(1, len(Rus))
for i in range(0, len(Rus)):
    for j in range(0, len(k_s)):
        #noramp_startup.generic_noramp_params["omega"]=omegas[i]
        #noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
        noramp_simulations.def_optim_list(["k_0", "Ru"])
        time_series1=noramp_simulations.test_vals([k_s[j], Rus[i]], "timeseries")
        voltages=noramp_simulations.define_voltages()
        #ax[i].plot(noramp_simulations.time_vec, voltages)
        plt.plot(voltages[val:], time_series1, label=str(k_s[j])+"_"+str(Rus[i]), linestyle=random_linestyles[i], color=random_colours[j])
    #ax[i].set_title('Ru='+ str(Rus[i]))
    #ax[i].set_xlabel("Voltage(V)")
    #ax[i].set_ylabel("Current(mA)")
    #ax[i].legend()
plt.legend()
plt.show()
