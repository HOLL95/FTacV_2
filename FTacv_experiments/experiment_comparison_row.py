import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
import math
from multiplotter import multiplot
from harmonics_plotter import harmonics
import time
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_6"
concs=["1e-1M", "1e0M"]
file_numbers=[str(x) for x in range(1, 4)]
#fig=multiplot(4, 4, **{"harmonic_position":3, "num_harmonics":5, "orientation":"portrait", "fourier_position":2, "plot_width":5})
#keys=sorted(fig.axes_dict.keys())
font = {'family' : 'normal',
        'size'   : 12}

plt.rc('font', **font)
letters=["A", "B", "C", "D", "E", "F", "G", "H"]
plot_counter=0
h_counter=0
pos=0
def RMSE(series1, series2):
    return np.sqrt((np.sum(1/(len(series1))*np.power(np.subtract(series1, series2),2))))

file="Noramp_2_cv_high_ru_alpha_disp"
method="timeseries"
noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
run="Run_2"
counter=1
ramped_file="Ramped_3_cv_high_ru.ts"
ramp_data_class=single_electron(("/").join([dir_path, results_dict,Electrode,"Run_2", ramped_file]), {}, {}, {}, {}, False)
#ramp_data_class.time_vec=ramp_data_class.other_values["experiment_time"]
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha_mean", "alpha_std"]
param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
param_vals[master_optim_list.index("Cdl")]=1e-5
ramp_param_vals=np.delete(param_vals, master_optim_list.index("cap_phase"))
ramp_param_vals[master_optim_list.index("CdlE1")]=0
ramp_param_vals[master_optim_list.index("CdlE2")]=0
ramp_param_vals[master_optim_list.index("Cdl")]=1e-5
noramp_results.def_optim_list(master_optim_list)
ramp_data_class.dim_dict["alpha_mean"]=0.5
ramp_data_class.dim_dict["alpha_std"]=1e-3
ramp_data_class.param_bounds=noramp_results.param_bounds
ramp_data_class.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","phase", "alpha_mean", "alpha_std"])
ramp_data_class.simulation_options["dispersion_bins"]=16
cmaes_time=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))*1e3
ramped_cmaes_time=ramp_data_class.i_nondim(ramp_data_class.test_vals(ramp_param_vals, method))*1e3
current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])#[0::dec_amount]
voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])#[0::dec_amount]
voltage_plots=noramp_results.e_nondim(noramp_results.define_voltages())
voltage_times=noramp_results.t_nondim(noramp_results.time_vec)
time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])#[0::dec_amount]
ramp_current_results=ramp_data_class.i_nondim(ramp_data_class.other_values["experiment_current"])#[0::dec_amount]
ramp_voltage_results=ramp_data_class.e_nondim(ramp_data_class.define_voltages())#[0::dec_amount]
ramp_time_results=ramp_data_class.t_nondim(ramp_data_class.time_vec)
num_harms=6
start_harm=1
end_harm=start_harm+num_harms
harms=harmonics(range(start_harm, end_harm),noramp_results.dim_dict["omega"] , 0.05)
noramp_harms=harms.generate_harmonics(time_results, cmaes_time)
ramp_harms=harms.generate_harmonics(ramp_time_results, ramped_cmaes_time)
all_harms=[ramp_harms, noramp_harms]
fourier_funcs=[np.real, np.imag]
fourier_ylabels=["Real", "Imaginary"]
fourier_times=[ramp_time_results, time_results]
fourier_currents=[ramped_cmaes_time, cmaes_time]
harm_xlabels=["Time(s)", "Voltage(V)"]
fig=multiplot(2, 4, **{"harmonic_position":3, "num_harmonics":num_harms, "fourier_position":2,"orientation":"portrait", "plot_width":6, "col_spacing":2})
keys=sorted(fig.axes_dict.keys())

fig.axes_dict["col1"][0].plot(ramp_time_results, ramp_voltage_results)
fig.axes_dict["col1"][0].set_xlabel("Time(s)")
fig.axes_dict["col1"][0].set_ylabel("Voltage(V)")
fig.axes_dict["col1"][1].plot(voltage_times, voltage_plots)
fig.axes_dict["col1"][1].set_xlabel("Time(s)")
fig.axes_dict["col1"][1].set_ylabel("Voltage(V)")
fig.axes_dict["col2"][0].plot(ramp_voltage_results, ramped_cmaes_time)
fig.axes_dict["col2"][0].set_xlabel("Time(s)")
fig.axes_dict["col2"][0].set_ylabel("Current(mA)")
fig.axes_dict["col2"][1].plot(voltage_results, cmaes_time)
fig.axes_dict["col2"][1].set_xlabel("Voltage(V)")
fig.axes_dict["col2"][1].set_ylabel("Current(mA)")
for i in range(0, 2):
    for j in range(0, 2):
        pos=(i*2)+j
        fig.axes_dict["col3"][pos].set_xlabel("Frequency(Hz)")
        fig.axes_dict["col3"][pos].set_ylabel(fourier_ylabels[j])
        if j==0:
            fig.axes_dict["col3"][pos].set_xticks([])
        harms.harmonic_selecter(fig.axes_dict["col3"][pos], fourier_currents[i],fourier_times[i], box=False, arg=fourier_funcs[j])
for i in range(0, 2):
    for j in range(0, num_harms):
        if i==0:
            x=ramp_time_results
        else:
            x=voltage_results
        pos=(i*num_harms)+j
        fig.axes_dict["col4"][pos].plot(x, all_harms[i][j,:]*1e3)
        twiny=fig.axes_dict["col4"][pos].twinx()
        twiny.set_ylabel(j+1, rotation=0)
        twiny.set_yticks([])
        if j==num_harms-1:
            fig.axes_dict["col4"][pos].set_xlabel(harm_xlabels[i])
        else:
            fig.axes_dict["col4"][pos].set_xticks([])
        if j==num_harms//2:
            fig.axes_dict["col4"][pos].set_ylabel("Current($\\mu A$)")
letter_count=0
y_pos=[1.1, 1.1, 1.195, 1.58]
for i in range(0, len(keys)):
    letter_idx=[0, len(fig.axes_dict[keys[i]])//2]

    for j in range(0, len(letter_idx)):
        pos=i+(j*(len(keys)))
        fig.axes_dict[keys[i]][letter_idx[j]].text(-0.1, y_pos[i], letters[pos], transform=fig.axes_dict[keys[i]][letter_idx[j]].transAxes,
          fontsize=14, fontweight='bold', va='top', ha='right')
        letter_count+=1

fig=plt.gcf()
fig.set_size_inches((14, 9))
#fig.tight_layout()
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.94, wspace=0.29, hspace=0.2)
plt.show()
save_path="experiment_comparison.png"
fig.savefig(save_path, dpi=500)
#plt.subplots_adjust(left =0.05 , right=0.97, wspace=0.20, top=0.91, bottom=0.06, hspace=0.05)
plt.show()
