import numpy as np
import matplotlib.pyplot as plt
import math
import time
start=time.time()
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
print "import", time.time()-start
types=["current", "voltage"]
start=time.time()
noramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={"GC4_1_cv":types, "GC4_2_cv":types, "GC4_3_cv":types}, dec_amount=4)
ramp_startup=FTACV_initialisation(experimental_fitting=False, file_dict={"GC4_1_ramp_cv":types,}, dec_amount=64)
print "read", time.time()-start
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":False,
    "dispersion":False,
    "dispersion_bins":20,
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
    "dispersion_bins":20,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"ramped",
    "label": "MCMC",
    "optim_list":[]
}
noramp_other_values={
    "filter_val": 0.5,
    "harmonic_range":range(3,7,1),
    "experiment_time": False,#noramp_startup.time_results["GC4_1_cv"],
    "experiment_current":False, #noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":False,#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":int(2e4),
}

ramped_other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,9,1),
    "experiment_time":False, #ramp_startup.time_results["GC4_1_ramp_cv"],
    "experiment_current":False, #ramp_startup.current_results["GC4_1_ramp_cv"],
    "experiment_voltage":False,#ramp_startup.voltage_results["GC4_1_ramp_cv"],
    "bounds_val":20,
    "signal_length":int(5e5),
}

noramp_simulations=single_electron(noramp_startup.generic_noramp_params, simulation_options, noramp_other_values)
ramped_simulations=single_electron(ramp_startup.generic_ramped_params, ramped_simulation_options, ramped_other_values)
harm_class=harmonics(noramp_other_values["harmonic_range"],noramp_simulations.nd_param.omega, 0.2)
ramp_ts=noramp_simulations.i_nondim(ramped_simulations.test_vals([], "timeseries"))*1e3
noramp_ts=noramp_simulations.i_nondim(noramp_simulations.test_vals([], "timeseries"))*1e3
volts_noramp=noramp_simulations.e_nondim(noramp_simulations.define_voltages())                   #noramp_simulations.other_values["experiment_voltage"])
volts_ramp=ramped_simulations.e_nondim(ramped_simulations.define_voltages())
volts_ramp=ramped_simulations.e_nondim(ramped_simulations.define_voltages())*1e3
volts_noramp=noramp_simulations.e_nondim(noramp_simulations.define_voltages())*1e3
time_noramp=noramp_simulations.t_nondim(noramp_simulations.time_vec)
time_ramp=ramped_simulations.t_nondim(ramped_simulations.time_vec)
noramp_harm=harm_class.generate_harmonics(time_noramp, noramp_ts)
ramp_harm=harm_class.generate_harmonics(time_ramp, ramp_ts)
harm_len=2
normal_plots=3
normal_axes=[]
harm_axes=[]
plt.rcParams.update({'font.size': 12})
xlabels=["Time(s)", "Time(s)", "Frequency(Hz)"]
ylabels=["Voltage(mV)", "Current(mA)", "Amplitude"]
letters=["A", "B", "C", "E", "F", "G"]
fig=plt.figure()
total_plots=4
cols=2
spacing=2
plot_height=2
for j in [0,1]:
    for i in range(0, normal_plots):
        normal_axes.append(plt.subplot2grid((1+spacing+normal_plots+harm_class.num_harmonics*total_plots*plot_height,cols), (i*((plot_height*harm_class.num_harmonics)+spacing),j), rowspan=harm_class.num_harmonics*plot_height, colspan=1))
        normal_axes[i].text(-0.1, 1.15, letters[i], transform=normal_axes[i].transAxes,
          fontsize=14, fontweight='bold', va='top', ha='right')
        normal_axes[i].set_xlabel(xlabels[i])
        normal_axes[i].set_ylabel(ylabels[i])
for i in range(0, normal_plots):
    normal_axes[i+normal_plots].set_xlabel(xlabels[i])
    normal_axes[i+normal_plots].set_ylabel(ylabels[i])
    normal_axes[i+normal_plots].text(-0.1, 1.15, letters[i+normal_plots], transform=normal_axes[i+normal_plots].transAxes,
      fontsize=14, fontweight='bold', va='top', ha='right')
counter=0
letter_counter=-1
harm_letters=["D", "H"]
for j in [0,1]:
    for i in range(0, harm_class.num_harmonics):
        print (1+harm_class.num_harmonics)*normal_plots+i*2
        harm_axes.append(plt.subplot2grid((1+spacing+normal_plots+harm_class.num_harmonics*total_plots*plot_height,cols), ((normal_plots*(plot_height*harm_class.num_harmonics+spacing))+(i*plot_height),j), rowspan=plot_height, colspan=1))

        if i==0:
            print "hey"
            letter_counter+=1
            harm_axes[counter].text(-0.1, 1.15, harm_letters[letter_counter], transform=harm_axes[counter].transAxes,
              fontsize=14, fontweight='bold', va='top', ha='right')
        counter+=1
#harm_axes.append(plt.subplot2grid((harm_class.num_harmonics,harm_len*2), (0,0), colspan=harm_len))
#harm_axes[i].yaxis.set_label_position("right")
#harm_axes[i].set_ylabel(str(harm_class.harmonics[i]), rotation=0)


normal_axes[0].plot(time_noramp, volts_noramp)
#normal_axes[0].set_title("Sinusoidal voltage vs time")
normal_axes[1].plot(time_noramp, noramp_ts)
#normal_axes[1].set_title("Sinusoidal current vs time")
harm_class.harmonic_selecter(normal_axes[2], noramp_ts, time_noramp)
#normal_axes[2].set_title("Sinusoidal absolute Fourier transform")
normal_axes[3].plot(time_ramp, volts_ramp)
#normal_axes[3].set_title("Ramped voltage vs time")
normal_axes[4].plot(time_ramp, ramp_ts)
#normal_axes[4].set_title("Ramped current vs time")
harm_class.harmonic_selecter(normal_axes[5], ramp_ts, time_ramp)
#normal_axes[5].set_title("Ramped absolute Fourier transform")
harm_counter=0
#harm_axes[harm_counter].set_title("Sinusoidal harmonics")
#normal_axes[6].text(1.17, 1.15, "D", transform=normal_axes[2].transAxes,
 # fontsize=14, fontweight='bold', va='top', ha='right')
#normal_axes[7].text(1.17, 1.15, "H", transform=normal_axes[5].transAxes,
#    fontsize=14, fontweight='bold', va='top', ha='right')

for i in range(0, len(noramp_harm)):
    harm_axes[harm_counter].plot(time_noramp, noramp_harm[i,:]*1e3)
    harm_axes[harm_counter].yaxis.set_label_position("right")
    harm_axes[harm_counter].set_ylabel(str(harm_class.harmonics[i]), rotation=0)
    if (i+1)==len(noramp_harm):
        harm_axes[harm_counter].set_xlabel("Time(s)")
    else:
        harm_axes[harm_counter].set_xticks([])

    harm_counter+=1
#harm_axes[harm_counter].set_title("Ramped harmonics")

for i in range(0, len(ramp_harm)):
    harm_axes[harm_counter].plot(time_ramp, abs(ramp_harm[i,:])*1e3)
    harm_axes[harm_counter].yaxis.set_label_position("right")
    harm_axes[harm_counter].set_ylabel(str(harm_class.harmonics[i]), rotation=0)

    locs=harm_axes[harm_counter].get_yticks()
    print locs
    new_locs=np.linspace(locs[0], locs[-1], 5),2
    new_locs=[round(x,2) for x in new_locs[0]]
    harm_axes[harm_counter].set_yticks(new_locs[1:3])
    if (i+1)==len(noramp_harm):
        harm_axes[harm_counter].set_xlabel("Time(s)")
    else:
        harm_axes[harm_counter].set_xticks([])
    harm_counter+=1
plt.subplots_adjust(left =0.11 , right=0.96, wspace=0.27, top=0.95, bottom=0.05)
fig.text(0.52, 0.15, 'Current($\mu$A)', ha='center', va='center', rotation='vertical')
fig.text(0.03, 0.15, 'Current($\mu$A)', ha='center', va='center', rotation='vertical')

plt.show()
