import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import sys
import matplotlib.ticker as ticker
import math
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
                print(name)
                return np.loadtxt(root+"/"+Electrode+"/"+name, skiprows=skiprows)
def change_param(params, optim_list, parameter, value):
    param_list=copy.deepcopy(params)
    param_list[optim_list.index(parameter)]=value
    return param_list
def chain_appender(chains, param):
    new_chain=chains[0, :, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, :, param])
    return new_chain
def plot_params(titles, set_chain):
    for i in range(0, len(titles)):
        print(i)
        axes=plt.subplot(3,5,i+1)
        plot_chain=chain_appender(set_chain, i)
        if abs(np.mean(plot_chain))<0.001:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3e'))
        else:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        axes.hist(plot_chain, alpha=0.4,bins=20, stacked=True, edgecolor='black')
        lb, ub = axes.get_xlim( )
        axes.set_xticks(np.linspace(lb, ub, 4))
        axes.set_xlabel(titles[i])
        axes.set_ylabel('frequency')
        axes.set_title(graph_titles[i])
files=os.listdir('.')
dcv_results=("_").join([Electrode, "Electrode", "Direct_CV", str(3), str(3)])
dcv_blank=("_").join(["Blank", Electrode, "Electrode", "Direct_CV", str(3), str(3)])
exp_results=find(dcv_results, one_above, Electrode,1)
blank_results=find(dcv_blank, one_above, Electrode,1)
dcv_currents=(exp_results[:,2]-blank_results[:,2])
dcv_voltages=exp_results[:,1]
dcv_times=exp_results[:,0]
ramped_current=("_").join([Electrode, "Electrode", "Ramped","3", "cv", "current" ])
ramped_voltage=("_").join([Electrode, "Electrode", "Ramped","3", "cv", "voltage" ])
voltage_results=find(ramped_voltage, one_above, Electrode)
current_results=find(ramped_current, one_above, Electrode)
noramp_current=find("Yellow_Electrode_noramp_3_cv_current", one_above, Electrode)
noramp_voltage=find("Yellow_Electrode_noramp_3_cv_voltage", one_above, Electrode)
dec_amount=32
noramp_current_results=noramp_current[0::dec_amount, 1]
noramp_time_results=noramp_current[0::dec_amount, 0]
noramp_voltage_results=noramp_voltage[0::dec_amount, 1]
ramped_time=current_results[:,0]
ramped_current=current_results[:,1]
ramped_voltage=voltage_results[:,1]
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    "E_0":0.2,
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-6, #(capacitance parameters)
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
    'num_peaks': 50
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
time_start=2/(param_list["omega"])
simulation_options={
    "no_transient":time_start,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":30,
    "test": False,
    "method": "sinusoidal",
    "phase_only": False,
    "likelihood":likelihood_options[0],
    "numerical_method": solver_list[1],
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(3,9,1)),
    "experiment_time": noramp_time_results,
    "experiment_current": noramp_current_results,
    "experiment_voltage":noramp_voltage_results,
    "bounds_val":20,
    "signal_length":int(1e4)
}
param_bounds={
    'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.95*param_list['omega'],1.05*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [40, 1e4],  #     (uncompensated resistance ohms)
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
noramp_class=single_electron(None, param_list, simulation_options, other_values, param_bounds)
unit_dict={
    "E_0": "V",
    'E_start': "V", #(starting dc voltage - V)
    'E_reverse': "V",
    'omega':"Hz",#8.88480830076,  #    (frequency Hz)
    'd_E': "V",   #(ac voltage amplitude - V) freq_range[j],#
    'v': '$s^{-1}$',   #       (scan rate s^-1)
    'area': '$cm^{2}$', #(electrode surface area cm^2)
    'Ru': "$\\Omega$",  #     (uncompensated resistance ohms)
    'Cdl': "F", #(capacitance parameters)
    'CdlE1': "",#0.000653657774506,
    'CdlE2': "",#0.000245772700637,
    'CdlE3': "",#1.10053945995e-06,
    'gamma': 'mol cm^{-2}$',
    'k_0': 's^{-1}$', #(reaction rate s-1)
    'alpha': "",
    "E0_mean":"V",
    "E0_std": "V",
    "k0_shape":"",
    "k0_loc":"",
    "k0_scale":"",
    "cap_phase":"rads",
    'phase' : "rads",
    "":"",
    "noise":"",
}
fancy_names={
    "E_0": '$E^0$',
    'E_start': '$E_{start}$', #(starting dc voltage - V)
    'E_reverse': '$E_{reverse}$',
    'omega':'$\\omega$',#8.88480830076,  #    (frequency Hz)
    'd_E': "$\\Delta E$",   #(ac voltage amplitude - V) freq_range[j],#
    'v': "v",   #       (scan rate s^-1)
    'area': "Area", #(electrode surface area cm^2)
    'Ru': "Ru",  #     (uncompensated resistance ohms)
    'Cdl': "$C_{dl}$", #(capacitance parameters)
    'CdlE1': "$C_{dlE1}$",#0.000653657774506,
    'CdlE2': "$C_{dlE2}$",#0.000245772700637,
    'CdlE3': "$C_{dlE3}$",#1.10053945995e-06,
    'gamma': '$\\Gamma',
    'k_0': '$k_0', #(reaction rate s-1)
    'alpha': "$\\alpha$",
    "E0_mean":"$E^0 \\mu$",
    "E0_std": "$E^0 \\sigma$",
    "cap_phase":"Capacitance phase",
    'phase' : "Phase",
    "":"Experiment",
    "noise":"$\sigma$",
}
Titles={
    'omega':'Input frequency',#8.88480830076,  #    (frequency Hz)
    'd_E': "Amplitude",   #(ac voltage amplitude - V) freq_range[j],#
    'v': "Scan rate",   #       (scan rate s^-1)
    'area': "Area", #(electrode surface area cm^2)
    'Ru': "Uncompensated resistance",  #     (uncompensated resistance ohms)
    'Cdl': "Linear capacitance", #(capacitance parameters)
    'CdlE1': "First order capacitance",#0.000653657774506,
    'CdlE2': "Second order capacitance",#0.000245772700637,
    'CdlE3': "Third order capacitance",#1.10053945995e-06,
    'gamma': 'Surface coverage',
    'k_0': 'Rate constant', #(reaction rate s-1)
    'alpha': "Symmetry factor",
    "E0_mean":"Themodynamic mean",
    "E0_std": "Thermodynamic standard deviation",
    "cap_phase":"Capacitance phase",
    'phase' : "Phase",
    "":"Experiment",
    "noise":"Noise",
}
#f=open(filename, "r")
optim_list=['E0_mean',"E0_std",'k_0',"Ru","Cdl","CdlE1", "CdlE2",'gamma', "cap_phase","phase", "noise"]
titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
graph_titles=[Titles[x] for x in optim_list]
true_params=[-0.4, 1e1, 0.000134*1e4,20.0]
noise_vals=[0.01, 0.02, 0.005, 0]
counter=-1
extension=".txt"
chain_select=1
j=0
#fig=plt.figure(num=None, figsize=(12,9), dpi=120, facecolor='w', edgecolor='k')
filename_list=[]
number_list=[]
num_params=len(optim_list)
electrode="Yellow"
folder="MCMC_runs"
for i in range(1, 9):
    filename="Noramp_3_timeseriesMCMC_"+str(i)  +"_run1"
    #Run4 - Constrained K0 and Ru, Run three fixed Ru, Alpha (not timeseries) has constrained Ru and unconstrained k, Run 6 has fixed ru and unconstrained k, run 5 has unconstrained both
    chains=np.load(("/").join([electrode, folder, filename]))
#chains2=np.load('GC4_MCMC_1_low_ru')
    pints.plot.trace(chains)
    #
    #plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92, wspace=0.30, hspace=0.33)
#plot_params(titles, chains2)
    plt.show()
    #plot_params(titles,chains[:, 30000:, :])
    #plt.show()
    params=np.zeros(len(optim_list)-1)
    for i in range(0, len(optim_list)-1):
        params[i]=np.mean(chains[0, 30000:, i])
    noramp_class.def_optim_list(optim_list[:-1])
    print(list(params))
    test=noramp_class.test_vals(params, "timeseries")
    plt.plot(noramp_class.other_values["experiment_voltage"],test)
    plt.plot(noramp_class.other_values["experiment_voltage"],noramp_class.other_values["experiment_current"], alpha=0.7)

    plt.show()
