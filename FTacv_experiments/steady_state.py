import numpy as np
import matplotlib.pyplot as plt
import math
import time
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
import os
import pickle
import pints
dir_path = os.path.dirname(os.path.realpath(__file__))
types=["current", "voltage"]
exp="Experimental-221119"
bla="Blank-110919"
resistances=["high_ru", "low_ru", "fixed_ru"]
ru_upper_bound=[700, 1e4, 50]
ru_pick=0
resistance_type=resistances[ru_pick]
print(resistance_type)
exp_type=exp
if exp_type==bla:
    extra="Blank/"
else:
    extra=""
data_path="/experiment_data_2/"+exp_type
Electrode="Red"
folder="Noramp"
for lcv_1 in range(1, 3):
    Method =str(lcv_1)+"_cv"
    type="current"
    type2="voltage"
    path=("/").join([dir_path, data_path, folder, Electrode])
    files= os.listdir(path)
    for data in files:
        if (Method in data)  and (type in data):
            results=np.loadtxt(path+"/"+data)
            print(data)
        elif (Method in data)  and (type2 in data):
            voltages=np.loadtxt(path+"/"+data)

    dec_amount=32
    de=300e-3
    estart=260e-3-de
    ereverse=estart+2*de
    current_results1=results[0::dec_amount, 1]
    time_results1=results[0::dec_amount, 0]
    voltage_results1=voltages[0::dec_amount, 1]
    results_dict={"experiment_voltage": voltage_results1,
                    "experiment_time": time_results1,
                    "experiment_current": current_results1,}
    param_list={
        "E_0":0.2,
        'E_start': estart, #(starting dc voltage - V)
        'E_reverse': ereverse,
        'omega':8.94,#8.88480830076,  #    (frequency Hz)
        "original_omega":8.94,
        'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
        'v': 10.36e-3,   #       (scan rate s^-1)
        'area': 0.07, #(electrode surface area cm^2)
        'Ru': 1.0,  #     (uncompensated resistance ohms)
        'Cdl': 1e-5, #(capacitance parameters)
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
        'sampling_freq' : (1.0/200),
        'phase' : 3*(math.pi/2),
        "time_end": None,
        'num_peaks': int(15)
    }
    solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
    likelihood_options=["timeseries", "fourier"]
    time_start=1/(param_list["omega"])
    simulation_options={
        "no_transient":time_start,
        "numerical_debugging": False,
        "experimental_fitting":True,
        "dispersion":False,
        "dispersion_bins":16,
        "test": False,
        "method": "sinusoidal",
        "phase_only":False,
        "likelihood":likelihood_options[0],
        "numerical_method": solver_list[1],
        "label": "MCMC",
        "optim_list":[]
    }
    other_values={
        "filter_val": 0.5,
        "harmonic_range":list(range(3,9,1)),
        "experiment_time": time_results1,
        "experiment_current": current_results1,
        "experiment_voltage":voltage_results1,
        "bounds_val":20,
        "signal_length":int(1e4)
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
        'alpha': [0.3, 0.7],
        "cap_phase":[0, 2*math.pi],
        "E0_mean":[0.19, 0.25],
        "E0_std": [0.001, 0.2],
        "k0_shape":[0,2],
        "k0_loc":[0, 1e3],
        "k0_scale":[0,1e3],
        "k0_range":[1e2, 1e4],
        'phase' : [0, 2*math.pi]
    }
    noramp_fit=single_electron(None, param_list, simulation_options, other_values, param_bounds)
    noramp_fit.define_boundaries(param_bounds)
    time_results=noramp_fit.other_values["experiment_time"]
    current_results=noramp_fit.other_values["experiment_current"]
    voltage_results=noramp_fit.other_values["experiment_voltage"]
    peak_range=range(2, param_list["num_peaks"]+1)
    partition_idx=[tuple(np.where((time_results>=(x-1)) & (time_results<x))) for x in peak_range]
    ratios=np.zeros(len(peak_range))
    positions=np.zeros(len(peak_range))
    for j in range(1, len(partition_idx)):
        current_partition=current_results[partition_idx[j]]
        time_partition=time_results[partition_idx[j]]
        max_current_peak=max(current_partition)
        prev_max=max(current_results[partition_idx[j-1]])
        ratios[j-1]=max_current_peak/prev_max
        #print(time_results[np.where(time_results==max(current_results[partition_idx[j]]))])
        positions[j-1]=time_partition[np.where(current_partition==max_current_peak)]
    #print(positions)
    plt.plot(time_results, current_results, label="Data")
    plt.plot(time_results[0], current_results[0], color="Red", label="Peak amplitude ratio")
    plt.xlabel("Nondim time")
    plt.ylabel("Nondim current")
    plt.legend()
    ax=plt.gca()
    ax2=ax.twinx()
    ax2.plot(positions[:-1], ratios[:-1], color="red", marker="x")
    print(np.log(positions))
    coeffs=np.polyfit(np.log(positions[:-1]), ratios[:-1], 1)
    print(coeffs)
    regression=np.add(np.multiply(np.log(time_results), coeffs[0]),coeffs[1])
    #regression=np.add(regression, np.multiply(np.power(np.log(time_results),2), coeffs[0]) )
    #ax2.plot(time_results, regression, color="black", linestyle="--")
    ax2.axhline(1, color="black", linestyle="--")
    plt.show()
