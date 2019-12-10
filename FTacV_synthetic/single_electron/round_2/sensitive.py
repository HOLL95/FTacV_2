import numpy as np
import matplotlib.pyplot as plt
import math
import time
import sys
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
dirs=(dir_path.split("/"))
desired_folder="FTacv_experiments"
global_folder="FTacV_2"
desired_path=("/").join(np.append(dirs[1:dirs.index(global_folder)+1], desired_folder))
sys.path.append("/"+desired_path+"/")
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
import isolver_martin_brent
from harmonics_plotter import harmonics
import os
import pickle
import pints
import pints.plot
import sympy as sym
import copy
from SALib.analyze import sobol
dec_amount=32
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
sf_list=[200, 400, 800, 1000, 2000, 5000, 10000, 20000, 50000]
prev_data="data"
counter=0
sensitive_params=["E_0","Ru", "k_0", "alpha"]
#fig, axes=plt.subplots(len(sensitive_params), 2)
omega_range=np.arange(1,801 , 1)
k0_range=[1,50, 100,500,1000, 5000]#, 1000, 10000]
ru_range=[300]
sf_range=[200, 400, 1000, 2000, 4000, 8000]
for lcv_2 in sf_range:
    fisher_vals=np.zeros((len(sensitive_params), len(omega_range)))
    for lcv_1 in range(0, len(omega_range)):
        sf=lcv_2
        omega=omega_range[lcv_1]
        print(omega)
        param_list={
        "E_0":0.2,
        'E_start': estart, #(starting dc voltage - V)
        'E_reverse': ereverse,
        'omega':omega,#8.88480830076,  #    (frequency Hz)
        'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
        'v': 1e-3,   #       (scan rate s^-1)
        'area': 0.07, #(electrode surface area cm^2)
        'Ru': 1.0,  #     (uncompensated resistance ohms)
        'Cdl': 1e-6, #(capacitance parameters)
        'CdlE1': 0,#0.000653657774506,
        'CdlE2': 0,#0.000245772700637,
        'CdlE3': 0,#1.10053945995e-06,
        'gamma': 1e-10,
        "original_gamma":1e-10,        # (surface coverage per unit area)
        'k_0': lcv_2, #(reaction rate s-1)
        'alpha': 0.5,
        "E0_mean":0.2,
        "E0_std": 0.09,
        "k0_shape":0.954,
        "k0_loc":100,
        "k0_scale":50,
        "k0_range":1e3,
        "cap_phase":0,
        'sampling_freq' : (1.0/sf),
        'phase' : 3*(math.pi/2),
        "time_end": None,
        'num_peaks': 10
        }
        solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
        likelihood_options=["timeseries", "fourier"]
        time_start=5/(param_list["omega"])
        #start_idx=np.where(time_results1>time_start)
        simulation_options={
        "no_transient":False,
        "numerical_debugging": False,
        "experimental_fitting":False,
        "dispersion":False,
        "dispersion_bins":16,
        "test": False,
        "method": "sinusoidal",
        "likelihood":likelihood_options[0],
        "numerical_method": solver_list[1],
        "label": "MCMC",
        "phase_only":False,
        "optim_list":[]
        }
        other_values={
        "filter_val": 0.5,
        "harmonic_range":list(range(3,9,1)),
        "experiment_time": None,
        "experiment_current": None,
        "experiment_voltage":None,
        "bounds_val":20,
        "signal_length":int(2e4)
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
        'k_0': [10, 1e3], #(reaction rate lcv_2s-1)
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
        colours=["red", "green"]

        vs=[1.0/omega, 2*math.pi/omega, omega]
        T=(273+25)
        F=96485.3328959
        R=8.314459848
        E0=(R*T)/F
        #print(E0, "E0")
        #print(2.5/omega, omega)
        param_list["v"]=E0*omega
        syn_fit=single_electron(None, param_list, simulation_options, other_values, param_bounds)
        true_vals=[0.25, 100, 300, 1e-10, 1e-5, 0,0, 3*math.pi/2, 3*math.pi/2, 0.5]
        syn_fit.def_optim_list(["E_0", "k_0", "Ru", "gamma", "Cdl", "CdlE1", "CdlE2", "phase", "cap_phase", "alpha"])
        data2=syn_fit.test_vals(true_vals, "timeseries")
        voltages2=syn_fit.define_voltages()
        times2=syn_fit.time_vec
        after_transient=5/param_list["omega"]
        time_idx2=tuple(range(0, len(data2)))#tuple(np.where(times2>after_transient))
        #plt.subplot(1,2,1)

        plot_voltages=(voltages2)
        plot_current=(data2)
        #time_idx=tuple(range(0, len(voltages2)))
        time_idx=tuple(np.where(syn_fit.t_nondim(syn_fit.time_vec)>(1/omega)))
        """
        plt.subplot(1,2,1)
        plt.plot(syn_fit.time_vec, voltages2)
        plt.subplot(1,2,2)
        plt.title(omega)
        plt.plot(syn_fit.e_nondim(plot_voltages[time_idx]), syn_fit.i_nondim(plot_current[time_idx]), label=("v "+str(sf)))
        plt.show()
        """
        #plot_times=syn_fit.t_nondim(syn_fit.time_vec)
        variables=["Cdl", "dE", "Ru", "d_I", "gamma", "k_0", "theta", "alpha", "E", "I","E_0"]

        for j in [1e-3]:
            numerical_mat=np.zeros((len(times2), len(sensitive_params)))
            delta_val=j
            orig_params=copy.deepcopy(true_vals)

            for i in range(0, len(sensitive_params)):
                pos=syn_fit.optim_list.index(sensitive_params[i])
                orig_vals=copy.deepcopy(true_vals[pos])
                true_vals[pos]=orig_vals+(orig_vals*delta_val)
                #print(true_vals, "change")
                delta_ts_1=syn_fit.i_nondim(syn_fit.test_vals(true_vals, "timeseries"))
                #true_vals[pos]=orig_vals-(orig_vals*delta_val)
                delta_ts_2=syn_fit.i_nondim(plot_current)#syn_fit.test_vals(orig_params, "timeseries")
                #print(true_vals, "change2")
                derivative_val=np.divide((np.subtract(delta_ts_1, delta_ts_2)),(orig_vals*delta_val))
                numerical_mat[:, i]=derivative_val
                true_vals=copy.deepcopy(orig_params)
                #plt.subplot(len(sensitive_params), 1, i+1)
                #plt.title(sensitive_params[i])
                #plt.plot(derivative_val[time_idx], label="$\Delta$="+str(delta_val))
                #plt.legend()
        #plt.show()
        #plot_current=d_I
        #plt.plot(syn_fit.time_vec, abs(numerical_mat[:, i]), label=omega)
        #plt.subplot(1,2,1)
        #plt.plot(syn_fit.e_nondim(voltages2[time_idx]),plot_current[time_idx])#
        #plt.subplot(1,2,2)
        #plt.subplot(1,2,1)
        #plt.title("Ru")
        #plt.plot(syn_fit.time_vec, numerical_mat[:,1], label=str(lcv_3)+" "+str(lcv_2))#np.power(numerical_mat[:,1],2), label=str(omega)+" "+str(np.trapz(y=np.power(numerical_mat[:,1],2), x=(syn_fit.time_vec))))#
        #plt.legend()
        #plt.subplot(1,2,2)
        #plt.title("K_0")
        #plt.plot(syn_fit.time_vec, numerical_mat[:,2], label=str(lcv_2)+" "+str(lcv_3))

        for i in range(0, len(sensitive_params)):
            fisher_vals[i, lcv_1]=np.trapz(y=(np.power(numerical_mat[:, i][time_idx],2)))#, x=syn_fit.t_nondim(syn_fit.time_vec[time_idx]))

    """    for i in range(0, len(sensitive_params)):
            print(i)
            print(len(sensitive_params))

            ax=axes[i, counter]
            ax.set_xlabel("Time(s)")
            ax.set_ylabel("Current(A)")

            ax.plot(voltages2[time_idx], (plot_current[time_idx]), label="f="+str(omega))
            if i==0:
                ax.legend()
            ax2=ax.twinx()
            ax2.set_ylabel("dI/d"+"("+str(sensitive_params[i])+")")
            ax.set_title(sensitive_params[i])
            ax2.plot(voltages2[time_idx], ((numerical_mat[:, i][time_idx])), color=colours[counter], linestyle="--")
            print(sensitive_params[i], np.trapz(y=abs((numerical_mat[:, i][time_idx])), x=voltages2[time_idx]))
        counter+=1"""



    for i in range(0, len(sensitive_params)):
        plt.subplot(2, 2, i+1)
        plt.semilogy(omega_range, (fisher_vals[i, :]), label="$S. rate=$"+str(lcv_2))
        plt.xlabel("Omega(Hz)")
        plt.ylabel("Fisher information")
        plt.title("$"+sensitive_params[i]+"$")
        plt.legend()
plt.show()
