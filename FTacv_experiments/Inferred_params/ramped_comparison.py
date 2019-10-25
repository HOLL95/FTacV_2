import os
import sys
import numpy as np
import matplotlib.pyplot as plt
dir_path = os.path.dirname(os.path.realpath(__file__))
slash_idx=[i for i in range(len(dir_path)) if dir_path[i]=="/"]
one_above=dir_path[:slash_idx[-1]]
sys.path.insert(1, one_above)
from single_e_class_unified import single_electron
from harmonics_plotter import harmonics
Electrode="Yellow"
path=("/").join([dir_path , Electrode])
files=os.listdir(path)#
file_numbers=[]
counter=1
cols=["low", "high", "fixed"]
rows=["2"]
file1="Noramp_1_cv_low_ru.pkl"
ramped_current=("_").join([Electrode, "Electrode", "Ramped",rows[0], "cv", "current" ])
ramped_voltage=("_").join([Electrode, "Electrode", "Ramped",rows[0], "cv", "voltage" ])
def find(name, path, Electrode):
    for root, dirs, files in os.walk(path):
        if Electrode in dirs:
            files=os.listdir(root+"/"+Electrode)
            if name in files:
                print(name)
                return np.loadtxt(root+"/"+Electrode+"/"+name)
param_dict={}
for lcv_1 in range(0, len(cols)):
    for lcv_2 in range(0, len(rows)):
        filename=("_").join(["Noramp", rows[lcv_2], "cv", cols[lcv_1], "ru"])+"."
        print(filename)
        result=single_electron(path+"/"+filename)
        param_dict[rows[lcv_2]+"_"+cols[lcv_1]]=result.save_dict["params"][0]
voltage_results=find(ramped_voltage, one_above, Electrode)
current_results=find(ramped_current, one_above, Electrode)
time_results=current_results[:,0]
current_results=current_results[:,1]
ramped_params={
    "E_0":0.25,
    'E_start': -180e-3, #(starting dc voltage - V)
    'E_reverse': 620e-3,    #  (reverse dc voltage - V)
    'omega':8.94,  #    (frequency Hz)
    'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 29.8e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-5, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,
    "original_gamma":1e-10,          # (surface coverage per unit area)
    'k_0': 1.0, #(reaction rate s-1)
    "E0_mean":0.2,
    "E0_std": 0.09,
    "k0_shape":0.954,
    "k0_loc":100,
    "k0_scale":50,
    "k0_range":1e3,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    "cap_phase":0,
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 600
}

simulation_options_ramped={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":20,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"ramped",
    "label": "MCMC",
    "optim_list":[]
}

ramp_other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(2,6,1)),
    "experiment_time": time_results,
    "experiment_current":current_results, #noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":voltage_results[:,1],#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":len(current_results),
}
ramped_results=single_electron(None, ramped_params, simulation_options_ramped, ramp_other_values, result.param_bounds)
ramped_results.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"])
time_results=ramped_results.other_values["experiment_time"]
current_results=ramped_results.other_values["experiment_current"]
harm_class=harmonics(ramped_results.other_values["harmonic_range"], ramped_results.nd_param.omega*ramped_results.nd_param.c_T0, 0.05)
#results_dict={"1_experimental": current_results}
results_dict={}
for keys in list(param_dict.keys()):
    if "fixed" in keys:
        ramped_results.def_optim_list(["E0_mean", "E0_std","k_0","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"])
        ramped_results.dim_dict["Ru"]=65
    time_series=ramped_results.test_vals(param_dict[keys], "timeseries")
    results_dict[keys]=time_series
    ramped_results.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"])

data_harmonics=harm_class.generate_harmonics(time_results, current_results)
data_likelihood=harm_class.inv_objective_fun(ramped_results.kaiser_filter, current_results)
plot_dict_0={
    "experimental_harmonics":ramped_results.i_nondim(data_harmonics),
    "experimental_time_series":ramped_results.i_nondim(current_results),
}
print(ramped_results.dim_dict["phase"])
for keys in list(results_dict.keys()):

    current_time=results_dict[keys]
    print(len(current_time))
    plot_dict_0[keys[keys.index("_")+1:]+"zresistace_harmonics"]=ramped_results.i_nondim(harm_class.generate_harmonics(time_results, current_time))
    plot_dict_0[keys[keys.index("_")+1:]+"zresistace_time_series"]=ramped_results.i_nondim(current_time)
harm_class.harmonics_and_time(ramped_results.t_nondim(time_results),"Yellow2", "abs", **plot_dict_0)
#harm_class.inv_objective_fun(ramped_results.kaiser_filter, current_time))
