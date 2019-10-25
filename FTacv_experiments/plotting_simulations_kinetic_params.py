import numpy as np
import matplotlib.pyplot as plt
import math
import time
start=time.time()
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
import matplotlib.gridspec as gridspec
print("import", time.time()-start)
types=["current", "voltage"]
start=time.time()
noramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={}, dec_amount=4)
ramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={}, dec_amount=64)
print("read", time.time()-start)
val=500*0
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
    "signal_length":int(2e5),
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
#plt.rcParams.update({'font.size': 12})
noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
#noramp_simulations=single_electron(ramp_startup.generic_ramped_params, ramped_simulation_options, ramped_other_values)
harm_class=harmonics(noramp_other_values["harmonic_range"],noramp_simulations.nd_param.omega, 0.2)
noramp_simulations.def_optim_list(["k_0", "omega"])
omegas=[9, 9*25]
k_s=[50, 100, 150, 300, 600]
Rus=[0, 50, 100, 150, 300]
Cdls=[1e-6,1e-5, 1e-4]
#
#k_s=np.linspace(0.1, 0.3, 5)
fig = plt.figure(figsize=(9,9))


colors=["green", "gold", "red", "gold"]
alignment=["bottom", "left", "right", "top"]
outer_val=0.005
inner_val=0.05
plot_axes=[]
for i in range(3):
    plot_axes.append([])
    #outer
    outergs = gridspec.GridSpec(1, 1)
    outergs.update(bottom=0.1,left=(i%3)*.33+outer_val,
                   top=0.9,  right=(1+i%3)*.33-outer_val,
                   wspace=0.35, hspace=0.35)
    outerax = fig.add_subplot(outergs[0])
    outerax.set_title("$C_{dl}=$"+str(Cdls[i]))
    outerax.tick_params(axis='both',which='both',bottom=0,left=0,
                        labelbottom=0, labelleft=0)

    [outerax.spines[x].set_color(colors[i]) for x in alignment]
    [outerax.spines[x].set_linewidth(1.5) for x in alignment]
    #outerax.set_facecolor(colors[i])
    #outerax.patch.set_alpha(0.3)

    #inner
    gs = gridspec.GridSpec(2, 2)
    print((i%3)*.3+0.08, (1+i%3)*.3-0.05)
    gs.update(bottom=0.15,left=(i%3)*.33+inner_val,
                   top=0.85,  right=(1+i%3)*.33-(inner_val-outer_val),
                   wspace=0.35, hspace=0.35)
    for k in range(4):
        ax = fig.add_subplot(gs[k])
        plot_axes[i].append(ax)
        #ax.set_title('Axes Title {}'.format(k+1), color=colors[i])
#plt.show()
for q in range(0, len(Cdls)):
    for i in range(0, len(omegas)):
        k_axis=plot_axes[q][i]
        r_axis=plot_axes[q][i+2]
        for j in range(0, len(k_s)):
            noramp_startup.generic_noramp_params["omega"]=omegas[i]
            noramp_simulations=single_electron(None, noramp_startup.generic_noramp_params, simulation_options, noramp_other_values, param_bounds)
            first_cycles=1/omegas[i]
            noramp_simulations.time_idx=np.where(noramp_simulations.time_vec>first_cycles)
            noramp_simulations.time_idx=noramp_simulations.time_idx[0][0]
            noramp_simulations.def_optim_list(["k_0","Ru" ,"omega","Cdl"])
            time_series1=noramp_simulations.i_nondim(noramp_simulations.test_vals([k_s[j],0, omegas[i], Cdls[q]], "timeseries"))
            time_series2=noramp_simulations.i_nondim(noramp_simulations.test_vals([100,Rus[j], omegas[i], Cdls[q]], "timeseries"))
            voltages=noramp_simulations.define_voltages()
            #ax[i].plot(noramp_simulations.time_vec, voltages)
            voltages=noramp_simulations.e_nondim(voltages)
            k_axis.plot(voltages[noramp_simulations.time_idx:],time_series1*1e3, label="$k^0=$"+str(k_s[j]))#
            r_axis.plot(voltages[noramp_simulations.time_idx:],time_series2*1e3, label="$R_u=$"+str(Rus[j]))#
        k_axis.set_title('$\\omega$='+ str(omegas[i])+"$Hz$")
        k_axis.set_xlabel("Voltage(V)")
        k_axis.set_ylabel("Current(mA)")
        k_axis.legend()
        r_axis.set_title('$\\omega$='+ str(omegas[i])+"$Hz$")
        #ax[1,i].set_title('$\\omega$='+ str(omegas[i]/math.pi)+"$\pi$")
        r_axis.set_xlabel("Voltage(V)")
        r_axis.set_ylabel("Current(mA)")
        r_axis.legend()
plt.show()
