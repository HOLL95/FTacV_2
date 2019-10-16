import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy
import copy
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
experiment=str(1)
repeat=str(1)
file1="Noramp_1_cv_low_ru.pkl"


def find(name, path, Electrode, skiprows=0):
    for root, dirs, files in os.walk(path):
        if Electrode in dirs:
            files=os.listdir(root+"/"+Electrode)
            if name in files:
                print name
                return np.loadtxt(root+"/"+Electrode+"/"+name, skiprows=skiprows)
param_dict={}
filename_1=["Ramped_3_cv_high_ru.ts"]
param_results=[]
labels=["Ramped"]
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"]
optim_dict={}



for names in filename_1:
    result=single_electron(path+"/"+names)
    param_results.append([result.save_dict["params"][0][result.save_dict["optim_list"].index(key)] if  (key in result.save_dict["optim_list"]) else result.dim_dict[key] for key in master_optim_list])
exp_currents=[]
blank_currents=[]
for i in range(1, 4):
    for j in range(1, 4):
        #plt.subplot(2,3,i)
        dcv_results=("_").join([Electrode, "Electrode", "Direct_CV", str(i), str(j)])
        dcv_blank=("_").join(["Blank", Electrode, "Electrode", "Direct_CV", str(i), str(j)])
        exp_results=find(dcv_results, one_above, Electrode,1)
        blank_results=find(dcv_blank, one_above, Electrode,1)
        blank_currents.append(blank_results[:,2])
        exp_currents.append(exp_results[:,2]-blank_results[:,2])
voltages=exp_results[:,1]

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
    "experiment_time": exp_results[:,0],
    "experiment_current":exp_currents[4], #noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":exp_results[:,1],#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":len(exp_results[:,0]),
}
dcv_results=single_electron(None, dcv_params, simulation_options_dcv, dcv_other_values, result.param_bounds)
dcv_results.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"])
time_results=dcv_results.other_values["experiment_time"]
current_results=dcv_results.other_values["experiment_current"]
voltage_results=dcv_results.other_values["experiment_voltage"]
voltages2=dcv_results.define_voltages()

#f1=np.fft.fftfreq(len(current_results), time_results[1]-time_results[0])
#f2=np.fft.fftfreq(len(stored_current), stored_time[2]-stored_time[0])
#plt.plot(f1, np.fft.fft(current_results))
#plt.plot(f2, np.fft.fft(stored_current))
#plt.show()
#results_dict={"1_experimental": current_results}
results_dict={}
for i in range(0, len(param_results)):
    cdl_0=copy.deepcopy(param_results[i])
    for j in range(0, len(dcv_results.optim_list)):
        if "Cdl" in dcv_results.optim_list[j]:
            cdl_0[j]=0
    time_series=dcv_results.test_vals(cdl_0, "timeseries")
    print dcv_results.dim_dict["Cdl"]
    cdl_idx=tuple(np.where((voltages<0.16) | (voltages>0.38)))
    farad_idx=tuple(np.where((voltages>0.16) & (voltages<0.38)))
    cdl_only=copy.deepcopy(current_results)
    cdl_only[farad_idx]=0
    inter=np.interp(time_results[farad_idx], time_results, current_results)
    plt.plot(inter)
    plt.show()
    halfway=np.where(voltages==max(voltages))


    cdl_part_1=cdl_only[:halfway[0][0]]
    cdl_part_2=cdl_only[halfway[0][0]:]
    cs=scipy.interpolate.CubicSpline(half_time[0::32], cdl_part_1[0::32])
    cs_results=cs(half_time[0::32])
    cs_results=np.interp(half_time, half_time[0::32], cs_results)
    plt.plot(cs_results)
    plt.show()
    plt.plot(half_time[0::32], cdl_part_1-cs_results)
    plt.show()

    #cdl_only=np.interp(farad_idx, range(0, len(cdl_only)), cdl_only)
    plt.plot(cdl_part_1)
    plt.plot(cdl_part_2)
plt.show()
    #plt.plot(voltages[5:], time_series[5:])

#plt.plot(voltages,current_results)
#plt.show()

    #data_harms=harm_class.generate_harmonics(time_results, time_series)
    #exp_harms=harm_class.generate_harmonics(time_results, current_results)

#harm_class.plot_harmonics(time_results, "abs", exp=exp_harms,sim=data_harms)
#data_harmonics=harm_class.generate_harmonics(time_results, current_results)
#data_likelihood=harm_class.inv_objective_fun(dcv_results.kaiser_filter, current_results)
#results_dict["experimental"]=dcv_results.i_nondim(current_results)
#harm_class.harmonics_plus("Yellow2", "abs",time_results, **results_dict)
#harm_class.inv_objective_fun(dcv_results.kaiser_filter, current_time))
#"voltages": dcv_results["experiment_voltage"]
