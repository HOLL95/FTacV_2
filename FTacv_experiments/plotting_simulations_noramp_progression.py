import numpy as np
import matplotlib.pyplot as plt
import math
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
types=["current", "voltage"]
startup=FTACV_initialisation(experimental_fitting=True, file_dict={"GC4_1_cv":types, "GC4_2_cv":types, "GC4_3_cv":types}, dec_amount=4)
simulation_options={
    "no_transient":2000,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":40,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"sinusoidal",
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,9,1),
    "experiment_time": startup.time_results["GC4_1_cv"],
    "experiment_current": startup.current_results["GC4_1_cv"],
    "experiment_voltage":startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":int(2e4),
}

noramp_simulations=single_electron(startup.generic_noramp_params, simulation_options, other_values)
other_values["signal_length"]=int(5e5)
#ramped_simulations=single_electron(startup.generic_ramped_params, simulation_options, other_values)
noramp_simulations.def_optim_list(["E_0", "k_0","Ru","Cdl","CdlE1","CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])
no_disp=[0.2612884710818322, 77.8597602837287, 100.000000195054659, 3.4964907695564266e-05, 0.0906694301595769, -0.0030113656849529626, 1.3935113888379478e-10, 8.94157246298468, 4.926130603328263, 4.3535376944340864, 0.7896377851231787]
plot_no_disp=noramp_simulations.test_vals(no_disp, "timeseries")
noramp_simulations.def_optim_list(["E0_mean","E0_std", "k_0","Ru","Cdl","CdlE1","CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])
feasible_1=[0.23542237635148117, 0.004406518077359679, 1.6338863524591192, 80.0000000070073, 3.24102216543367e-05, 0.10683350182222259, -0.0041305496741190365, 1.531049485509304e-10, 8.94090414236465, 5.635236560276772, 4.327586365963918, 0.5420683576660721]
feasible_2=[0.2507124585192858, 0.012810161448688611, 81.25609725526283,100, 3.4455963089647964e-05, 0.08570032724521258, -0.0026762829980299324, 1.4591816196158548e-10, 8.940730783653084, 4.965641139765749, 4.375333040247454, 0.7811192097028208]
feasible_3=[0.2546334543150689, 0.03323126819215362, 116.24098132960108,144.78687238133318,3.4422425270083556e-05, 0.09197677289665342, -0.002971108200729257, 1.4130483196036997e-10,8.94090414236465, 4.947246577563367, 4.361083156674927, 0.5]
plot_feasible_1=noramp_simulations.test_vals(feasible_1, "timeseries")
plot_feasible_2=noramp_simulations.test_vals(feasible_2, "timeseries")
plot_feasible_3=noramp_simulations.test_vals(feasible_3, "timeseries")
noramp_simulations.def_optim_list(["E_0", "k_0","Ru","Cdl","CdlE1","CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])
no_cap_phase=[0.26785228085275314, 9999.999999251795, 1.000000000911962, 4.35578395627439e-05, 0.14999999999966418, -0.006762504610321346, 3.990440678447067e-11, 8.940830114295643, 4.537306794394063, 4.537306794394063, 0.8999999999708405]
plot_no_cap_phase=noramp_simulations.test_vals(no_cap_phase, "timeseries")
noramp_simulations.def_optim_list(["E_0", "k_0","Ru","Cdl",'gamma', 'omega',"phase","cap_phase","alpha"])
noramp_simulations.dim_dict["CdlE1"]=0
noramp_simulations.dim_dict["CdlE2"]=0
noramp_simulations.dim_dict["CdlE3"]=0
linear_capacitance=[0.2678506598727337, 854.1054279611175, 199.9999999997499, 4.9609877512277555e-05, 2.009415004494725e-10, 8.940937614608458, 4.8724923726244205, 4.379359117553546, 0.1000000000223103]
plot_linear_cap=noramp_simulations.test_vals(linear_capacitance, "timeseries")
noramp_simulations.def_optim_list(["E_0", "k_0","Ru","Cdl","CdlE1",'gamma', 'omega',"phase","cap_phase","alpha"])
cdle1_cap=[0.2635909320634533, 94.15018823936819, 199.99999999905086, 3.776905397954239e-05, 0.027260979368400273, 1.6724687162517318e-10, 8.942016236647628, 4.961900510168174, 4.379713058778126, 0.823948010049349]
plot_cdle1_cap=noramp_simulations.test_vals(cdle1_cap, "timeseries")
plot_array=[plot_no_cap_phase,plot_linear_cap, plot_cdle1_cap,plot_feasible_1, plot_feasible_2, plot_feasible_3, plot_feasible_2]
num_plots=len(plot_array)
positions=[6,7,8,4,9,14,10]
title_pos=[6,7,8,4,10]
titles=["Initial attempt", "Independent capacitance phase", "Higher order capacitance", "Dispersion", "Experimental information"]
pos=zip(title_pos, titles)
pos_idx=0
plt.rcParams.update({'font.size': 10})
volts=noramp_simulations.e_nondim(noramp_simulations.other_values["experiment_voltage"]*1e3)
currs=noramp_simulations.i_nondim(noramp_simulations.other_values["experiment_current"]*1e3)
letters=["A", "B", "C", "D", "E", "F", "G"]
for i in range(0, len(plot_array)):
    plt.subplot(3, 5, positions[i])
    plt.xlabel("Voltage(mV)")
    plt.ylabel("Current(mA)")
    plt.plot(volts,currs)
    plt.plot(volts, noramp_simulations.i_nondim(plot_array[i])*1e3)
    ax=plt.gca()
    ax.text(-0.1, 1.15, letters[i], transform=ax.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')
    if positions[i]==pos[pos_idx][0]:
        plt.title(pos[pos_idx][1])
        pos_idx+=1
plt.show()
