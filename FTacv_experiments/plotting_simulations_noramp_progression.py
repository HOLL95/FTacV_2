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

linear_capacitance=[0.2678506598727337, 154.1054279611175, 50.9999999997499, 4.9609877512277555e-05, 1.009415004494725e-10, 8.940959682968415, 3*math.pi/2, 3*math.pi/2, 0.5000000000223103]

noramp_simulations.def_optim_list(["E_0", "k_0","Ru","Cdl",'gamma', 'omega',"phase","cap_phase","alpha"])
no_cap_phase=[0.26785228085275314, 9999.999999251795, 1.000000000911962, 4.35578395627439e-05, 0.14999999999966418, -0.006762504610321346, 3.990440678447067e-11, 8.940830114295643, 4.537306794394063, 4.537306794394063, 0.8999999999708405]
plot_no_cap_phase=noramp_simulations.test_vals(linear_capacitance, "timeseries")
noramp_simulations.def_optim_list(["E_0", "k_0","Ru","Cdl","CdlE1","CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])
no_disp=[0.2612884710818322, 77.8597602837287, 100.000000195054659, 3.4964907695564266e-05, 0.0906694301595769, -0.0030113656849529626, 1.3935113888379478e-10, 8.94157246298468, 4.926130603328263, 4.3535376944340864, 0.7896377851231787]
plot_no_disp=noramp_simulations.test_vals(no_disp, "timeseries")
noramp_simulations.def_optim_list(["E_0", "k_0","Ru","Cdl","CdlE1","CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])

independent_phase=[0.2546334543150689, 116.24098132960108,144.78687238133318,3.4422425270083556e-05, 0.09997677289665342, -0.002971108200729257, 9.5130483196036997e-11,8.940959682968415, 4.947246577563367, 4.361083156674927, 0.5]
no_independent_phase=[0.27545349898901317, 499.9999997853943, 1.9888896435113728e-08, 4.4045733036758575e-05, 0.149999999980652, -0.006907122380785267, 4.1163615775576164e-11, 8.941581035878608, 4.526707361403333, 4.526707361403333,0.7999999999239]
all_together=[0.23953731576330606, 0.06251048824509087, 54.80881131862505, 1.0, 3.401227622853829e-05, 0.09658746281397658, -0.0035025931533504253, 1.4039829708596723e-10, 8.940959682968415, 4.88941265694472, 4.331036700670979, 0.5725681192157199]
disped=[0.2546334543150689, 0.06323126819215362, 116.24098132960108,144.78687238133318,3.4422425270083556e-05, 0.09197677289665342, -0.002971108200729257, 1.4130483196036997e-10, 8.940959682968415, 3*math.pi/2, 3*math.pi/2, 0.5]
#disped=[0.19500581798491548, 0.06412642660733478, 1.0000000000000577, 499.99999999836524, 1.97326127580925e-18, 0.07266916327583202, 0.005288084406834424, 9.148401675280261e-10, 8.941581035878608, 4.526707361403333, 4.526707361403333, 0.6790813967967773]

plot_cdle1_cap=noramp_simulations.test_vals(no_independent_phase, "timeseries")
plot_independent_phase=noramp_simulations.test_vals(independent_phase, "timeseries")
noramp_simulations.def_optim_list(["E0_mean","E0_std", "k_0","Ru","Cdl","CdlE1","CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])

plot_all_together=noramp_simulations.test_vals(all_together, "timeseries")

plot_just_disp=noramp_simulations.test_vals(disped, "timeseries")
plot_array=[plot_no_cap_phase, plot_cdle1_cap, plot_independent_phase, plot_just_disp, plot_all_together]
num_plots=len(plot_array)
positions=[6,7,8,4,9,14,10]
title_pos=[6,7,8,4,10]
titles=["Initial attempt", "Higher order capacitance", "Independent capacitance phase", "Thermodynamic dispersion","Final fit"]
pos=zip(title_pos, titles)
pos_idx=0
plt.rcParams.update({'font.size': 10})
volts=noramp_simulations.e_nondim(noramp_simulations.other_values["experiment_voltage"]*1e3)
currs=noramp_simulations.i_nondim(noramp_simulations.other_values["experiment_current"]*1e3)
letters=["A", "B", "C", "D", "E", "F", "G"]
for i in range(0, len(plot_array)):
    plt.subplot(1, 5, i+1)
    plt.xlabel("Voltage(mV)")
    plt.ylabel("Current(mA)")
    plt.title(titles[i])
    plt.plot(volts,currs)
    plt.plot(volts, noramp_simulations.i_nondim(plot_array[i])*1e3)
    ax=plt.gca()
    ax.text(-0.1, 1.15, letters[i], transform=ax.transAxes,
      fontsize=12, fontweight='bold', va='top', ha='right')
plt.show()
