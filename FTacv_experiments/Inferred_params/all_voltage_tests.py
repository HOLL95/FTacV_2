import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import math
import time
import math
import copy
dir_path = os.path.dirname(os.path.realpath(__file__))
slash_idx=[i for i in range(len(dir_path)) if dir_path[i]=="/"]
one_above=dir_path[:slash_idx[-1]]
sys.path.insert(1, one_above)
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
Electrode="Yellow"
path=("/").join([dir_path , Electrode])
files=os.listdir(path)#
def find(name, path, Electrode, skiprows=0):
    for root, dirs, files in os.walk(path):
        if Electrode in dirs:
            files=os.listdir(root+"/"+Electrode)
            if name in files:
                print name
                return np.loadtxt(root+"/"+Electrode+"/"+name, skiprows=skiprows)
def change_param(params, optim_list, parameter, value):
    param_list=copy.deepcopy(params)
    param_list[optim_list.index(parameter)]=value
    return param_list
dcv_results=("_").join([Electrode, "Electrode", "Direct_CV", str(2), str(2)])
dcv_blank=("_").join(["Blank", Electrode, "Electrode", "Direct_CV", str(2), str(2)])
exp_results=find(dcv_results, one_above, Electrode,1)
blank_results=find(dcv_blank, one_above, Electrode,1)
dcv_currents=(exp_results[:,2]-blank_results[:,2])
dcv_voltages=exp_results[:,1]
dcv_times=exp_results[:,0]
cdl_idx=tuple(np.where((dcv_voltages<0.05) | (dcv_voltages>0.48)))
cdl_times=dcv_times[cdl_idx]
cdl_current=dcv_currents[cdl_idx]
cdl_voltage=dcv_voltages[cdl_idx]
interped_background=np.interp(dcv_times, cdl_times, cdl_current)
dcv_data=np.subtract(dcv_currents,interped_background)
ramped_current=("_").join([Electrode, "Electrode", "Ramped","3", "cv", "current" ])
ramped_voltage=("_").join([Electrode, "Electrode", "Ramped","3", "cv", "voltage" ])
voltage_results=find(ramped_voltage, one_above, Electrode)
current_results=find(ramped_current, one_above, Electrode)
ramped_time=current_results[:,0]
ramped_current=current_results[:,1]
ramped_voltage=voltage_results[:,1]
dcv_params={
    "E_0":0.25,
    'E_start': -0.15, #(starting dc voltage - V)
    'E_reverse': 0.6,    #  (reverse dc voltage - V)
    'omega':8.94,  #    (frequency Hz)
    'd_E': 0,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 30e-3,   #       (scan rate s^-1)
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
simulation_options_dcv={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":20,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"dcv",
    "label": "MCMC",
    "optim_list":[]
}
dcv_other_values={
    "filter_val": 0.5,
    "harmonic_range":range(2,6,1),
    "experiment_time": dcv_times,
    "experiment_current":dcv_data, #noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":dcv_voltages,#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":len(dcv_voltages),
}

param_names=["Noramp_3_cv_fixed_ru.fixed_alpha", "Ramped_3_cv_high_ru.ts", "Noramp_3_cv_high_ru.fourier"]
ramped_class=single_electron(path+"/"+param_names[1])
noramp_ts_class=single_electron(path+"/"+param_names[0])
noramp_f_class=single_electron(path+"/"+param_names[2])
dcv_class=single_electron(None, dcv_params, simulation_options_dcv, dcv_other_values, noramp_ts_class.param_bounds)
noramp_current=noramp_ts_class.i_nondim(noramp_ts_class.other_values["experiment_current"])
noramp_voltage=noramp_ts_class.e_nondim(noramp_ts_class.other_values["experiment_voltage"])
noramp_time=noramp_ts_class.t_nondim(noramp_ts_class.other_values["experiment_time"])
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","cap_phase","phase","alpha"]
class_list=[noramp_ts_class, ramped_class, noramp_f_class]
simulation_classes=[noramp_ts_class, ramped_class, dcv_class]
Title_list=["Pinning", "Ramped", "Fourier"]
param_list=[]
for i in range(0, len(class_list)):
    result=class_list[i]
    param_list.append([result.save_dict["params"][0][result.save_dict["optim_list"].index(key)] if  (key in result.save_dict["optim_list"]) else result.dim_dict[key] for key in master_optim_list])
for i in range(0, len(param_list)):
    print param_list
for i in range(0, len(simulation_classes)):
    simulation_classes[i].def_optim_list(master_optim_list)

noramp_simulation_list=[]
ramped_simulation_list=[]
dcv_simulation_list=[]
noramp_ts_class.dim_dict["omega"]=8.940641321267664
ramped_class.dim_dict["omega"]=8.884904955014296
ramp_harm=harmonics(range(2, 6), ramped_class.nd_param.omega, 0.1)
noramp_harm=harmonics(range(2, 6), noramp_ts_class.nd_param.omega, 0.1)
for i in range(0, len(param_list)):
    dcv_params=change_param(copy.deepcopy(param_list[i]), master_optim_list, "Cdl", 0)
    if i==1:
        param_list[i]=change_param(param_list[i], master_optim_list, "cap_phase", 3*math.pi/2)
        param_list[i]=change_param(param_list[i], master_optim_list, "phase", 3*math.pi/2)

    noramp_simulation_list.append(noramp_ts_class.i_nondim(noramp_ts_class.test_vals(param_list[i], "timeseries")))
    ramped_simulation_list.append(ramped_class.i_nondim(ramped_class.test_vals(param_list[i], "timeseries")))
    dcv_simulation_list.append(dcv_class.i_nondim(dcv_class.test_vals(dcv_params, "timeseries")))


num_harmonics=4
plot_width=3
num_runs=3
num_plots=5
seperation=1
rows=num_plots*num_harmonics+(num_plots-1)
cols=plot_width*num_runs+(num_runs-1)
ts_row_idx=np.multiply(range(0, 3),((num_harmonics+seperation)*2))
harm_idx=np.multiply(range(1, 4, 2), (num_harmonics+seperation))
sinusoidal_axes=[]
ramped_axes=[]
dcv_axes=[]
sinusoidal_harm_axes=[]
ramped_harm_axes=[]

for i in range(0, num_runs):
    sinusoidal_axes.append(plt.subplot2grid((rows, cols), (ts_row_idx[0], i*(seperation+plot_width)), rowspan=num_harmonics, colspan=plot_width))
    ramped_axes.append(plt.subplot2grid((rows, cols), (ts_row_idx[1], i*(seperation+plot_width)), rowspan=num_harmonics, colspan=plot_width))
    dcv_axes.append(plt.subplot2grid((rows, cols), (ts_row_idx[2], i*(seperation+plot_width)), rowspan=num_harmonics, colspan=plot_width))
for i in range(0, num_runs):
    for j in range(0, num_harmonics):
        sinusoidal_harm_axes.append(plt.subplot2grid((rows, cols), (harm_idx[0]+j, i*(seperation+plot_width)), rowspan=1, colspan=plot_width))
        print harm_idx[0]+j, i*(seperation+plot_width)
        ramped_harm_axes.append(plt.subplot2grid((rows, cols), (harm_idx[1]+j, i*(seperation+plot_width)), rowspan=1, colspan=plot_width))
noramp_harmonics=[]
noramp_data_harms=noramp_harm.generate_harmonics(noramp_time, noramp_current)
ramped_harmonics=[]
ramped_data_harms=ramp_harm.generate_harmonics(ramped_time, ramped_current)
for i in range(0, num_runs):
    noramp_harmonics.append(noramp_harm.generate_harmonics(noramp_time, noramp_simulation_list[i]))
    ramped_harmonics.append(ramp_harm.generate_harmonics(ramped_time, ramped_simulation_list[i]))
    if i==1:
        ramped_simulation_list[i]=np.append(np.fft.ifft(ramped_class.kaiser_filter(ramped_simulation_list[i])), [0,0])
        #ramp_harm.inv_objective_fun(ramped_class.kaiser_filter, ramped_simulation_list[i])
        #noramp_simulation_list[i]=ramp_harm.inv_objective_fun(noramp_ts_class.kaiser_filter, noramp_simulation_list[i])


    sinusoidal_axes[i].plot(noramp_voltage, noramp_simulation_list[i])
    sinusoidal_axes[i].plot(noramp_voltage, noramp_current, alpha=0.5)
    ramped_axes[i].plot(ramped_time, ramped_simulation_list[i])
    if i==1:
        ramped_axes[i].plot(ramped_time, np.append(np.fft.ifft(ramped_class.kaiser_filter(ramped_current)), [0,0]), alpha=0.5)
    else:
        ramped_axes[i].plot(ramped_time, ramped_current, alpha=0.5)
    dcv_axes[i].plot(dcv_voltages[5:], dcv_simulation_list[i][5:])
    dcv_axes[i].plot(dcv_voltages, dcv_data, alpha=0.5)
counter=0
for i in range(0, len(noramp_harmonics)):
    for j in range(0, num_harmonics):
        sinusoidal_harm_axes[counter].plot(noramp_time, noramp_harmonics[i][j,:])
        sinusoidal_harm_axes[counter].plot(noramp_time, noramp_data_harms[j,:])
        ramped_harm_axes[counter].plot(ramped_time, abs(ramped_harmonics[i][j,:]))
        ramped_harm_axes[counter].plot(ramped_time, abs(ramped_data_harms[j,:]))
        counter+=1
plt.show()
results_dict={}
results_dict["experimental"]=ramped_current
results_dict["simulation"]=ramped_simulation_list[1]
ramp_harm.harmonics_plus("Yellow2", "abs",ramped_time, **results_dict)
for i in range(0,4):
    plt.subplot(2,2, i+1)
    plt.plot(abs(ramped_harmonics[1][i,:]))
    plt.plot(abs(ramped_data_harms[i,:]))
plt.show()
