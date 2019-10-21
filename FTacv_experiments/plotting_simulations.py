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
noramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={"GC4_1_cv":types, "GC4_2_cv":types, "GC4_3_cv":types}, dec_amount=4)
ramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={"GC4_1_ramp_cv":types,}, dec_amount=64)
print("read", time.time()-start)
simulation_options={
    "no_transient":False,
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
    "experiment_time": False,#noramp_startup.time_results["GC4_1_cv"],
    "experiment_current":False, #noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":False,#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":int(2e4),
}

ramped_other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(1,9,1)),
    "experiment_time":False, #ramp_startup.time_results["GC4_1_ramp_cv"],
    "experiment_current":False, #ramp_startup.current_results["GC4_1_ramp_cv"],
    "experiment_voltage":False,#ramp_startup.voltage_results["GC4_1_ramp_cv"],
    "bounds_val":20,
    "signal_length":int(5e5),
}

noramp_simulations=single_electron(noramp_startup.generic_noramp_params, simulation_options, noramp_other_values)
ramped_simulations=single_electron(ramp_startup.generic_ramped_params, ramped_simulation_options, ramped_other_values)
harm_class=harmonics(noramp_other_values["harmonic_range"],noramp_simulations.nd_param.omega, 0.2)
