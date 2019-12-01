import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
from harmonics_plotter import harmonics
import pints.plot
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_2"
counter=1
ramped_file="Ramped_3_cv_high_ru.ts"
ramped_data_class=single_electron(("/").join([dir_path, results_dict,Electrode,run, ramped_file]), {}, {}, {}, {}, False)
ramped_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
ramped_data_class.def_optim_list(ramped_optim_list)
ramped_data_class.dim_dict["v_nondim"]=True
ramped_data_class.dim_dict["CdlE3"]=0
ramped_data_class.dim_dict["phase"]=0
ramped_data_class.dim_dict["cap_phase"]=0
#ramped_data_class.dim_dict["omega"]=8.95
ramped_data_class.dim_dict["original_omega"]=8.94
ramped_data_class.simulation_options["dispersion_bins"]=200
ramped_data_path=("/").join([dir_path, "experiment_data_2", "Experimental-120919", "Ramped", "Yellow"])
start_harm=1
end_harm=6
ramp_data_harm_class=harmonics(list(range(start_harm, end_harm)), ramped_data_class.nd_param.omega, 0.1)
ramped_data_class.harmonic_range=list(range(0, 11))
num_exps=4
spacing=1
num_harmonics=end_harm-start_harm
total_rows=(num_exps*num_harmonics)+(num_exps-1)*spacing
col_width=3
num_plots=4
total_cols=(col_width*num_plots)+(num_plots-1)*spacing
time_axes=[]
sinusoidal_axes=[]
ramped_axes=[]
fft_axes=[]
for i in range(0, num_exps):
    time_axes.append(plt.subplot2grid((total_rows, total_cols), ((i*(num_harmonics+spacing)), 0), rowspan=num_harmonics, colspan=col_width))
for i in range(0, num_exps):
    for j in [0,int(np.ceil(num_harmonics/2))]:
        print((j+(i*(num_harmonics+spacing))))
        fft_axes.append(plt.subplot2grid((total_rows, total_cols), (j+(i*(num_harmonics+spacing)), col_width+spacing), rowspan=int(np.floor(num_harmonics/2)), colspan=col_width))

for i in range(0, num_exps):
    for j in range(0, num_harmonics):
        sinusoidal_axes.append(plt.subplot2grid((total_rows, total_cols), (j+(i*(spacing+num_harmonics)), 2*(col_width+spacing)), rowspan=1, colspan=col_width))
for i in range(0, num_exps):
    for j in range(0, num_harmonics):
        ramped_axes.append(plt.subplot2grid((total_rows, total_cols), (j+(i*(spacing+num_harmonics)), 3*(col_width+spacing)), rowspan=1, colspan=col_width))

def likelihood_func(harmonic_range, times, time_series, omega):
    filter_val=0.1
    frequency=np.fft.fftfreq(len(times), times[1]-times[0])
    hanning=np.hanning(len(times))
    time_series=np.multiply(time_series, hanning)
    desired_idx=np.where((frequency>((harmonic_range[0]*omega)-(omega*filter_val)))&(frequency<(harmonic_range[-1]*omega)+(omega*filter_val)))
    Y=(np.fft.fft(time_series))
    return Y[desired_idx],frequency[desired_idx]
method="timeseries"
counter=0
j=0
counter2=0
funcs=[np.real, np.imag]
labels=["Real", "Imaginary"]
for i in range(2, 2+num_exps):
    ramped_current_results=np.loadtxt(ramped_data_path+"/Yellow_Electrode_Ramped_"+str(i)+"_cv_current")
    ramped_voltage_results=np.loadtxt(ramped_data_path+"/Yellow_Electrode_Ramped_"+str(i)+"_cv_voltage")[:,1]
    ramped_time_results=ramped_current_results[:,0]
    ramped_current_results=ramped_current_results[:,1]
    ramped_data_harmonics=ramp_data_harm_class.generate_harmonics(ramped_time_results, ramped_current_results)
    file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,"Run_3", file]))
    noramp_harm=harmonics(list(range(start_harm, end_harm)), noramp_results.nd_param.omega, 0.1)
    param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    noramp_results.def_optim_list(master_optim_list)
    harmonic_range=list(range(2, 9))

    #print(param_vals)

    cmaes_time_series=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))



    current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])
    voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])
    time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])
    cmaes_times=time_results#noramp_results.t_nondim(noramp_results.time_vec[noramp_results.time_idx:])
    test_voltages=voltage_results
    filtered_exp, filtered_exp_freq=likelihood_func(harmonic_range, cmaes_times, cmaes_time_series, noramp_results.nd_param.omega)
    filtered_data, filtered_data_freq=likelihood_func(harmonic_range, time_results, current_results,  noramp_results.nd_param.omega)

    #test_voltages=np.interp(cmaes_times, time_results, voltage_results)
    ramped_param_vals=[noramp_results.dim_dict[x] for x in ramped_optim_list]
    ramped_param_vals[-1]=0.5
    #change_params=["omega", "phase", "cap_phase"]
    #change_vals=[8.95, 0, 0]
    #for q in range(0, len(change_vals)):
    #    ramped_param_vals[ramped_optim_list.index(change_params[q])]=change_vals[q]
    ramped_time_series=ramped_data_class.i_nondim(ramped_data_class.test_vals(ramped_param_vals, method))
    ramped_times=ramped_data_class.t_nondim(ramped_data_class.time_vec[ramped_data_class.time_idx:])
    ramped_harmonics=ramp_data_harm_class.generate_harmonics(ramped_times, ramped_time_series)
    sinusoid_data_harmonics=noramp_harm.generate_harmonics(time_results, current_results)
    sinusoid_harmonics=noramp_harm.generate_harmonics(cmaes_times, cmaes_time_series)
    for q in range(0, 2):
        fft_axes[counter2].plot(filtered_exp_freq, funcs[q](filtered_exp), label=labels[q]+" Simulation")
        fft_axes[counter2].plot(filtered_data_freq, funcs[q](filtered_data), label=labels[q]+" Data")
        fft_axes[counter2].set_xlabel("Frequency(Hz)")
        fft_axes[counter2].set_ylabel(labels[q])
        #fft_axes[counter2].legend()
        counter2+=1

    time_axes[counter].plot(test_voltages, cmaes_time_series, label="Simulation")
    time_axes[counter].plot(voltage_results, current_results, alpha=0.7, label="Data")
    time_axes[counter].legend()
    time_axes[counter].set_xlabel("Voltage(V)")
    time_axes[counter].set_ylabel("Current(A)")
    #time_axes[counter].set_title("Sinusoidal " +str(i))
    #plt.plot(ramped_times, ramped_data_class.define_voltages()[ramped_data_class.time_idx:])
    #plt.show()
    harm_counter=0
    for harm_counter in range(0, len(ramped_harmonics)):
        ramped_axes[j].plot(ramped_times, abs(ramped_harmonics[harm_counter,:]*1e6), label="Simulation")
        ramped_axes[j].plot(ramped_time_results, abs(ramped_data_harmonics[harm_counter,:]*1e6), label="Data")
        sinusoidal_axes[j].plot(cmaes_times, (sinusoid_harmonics[harm_counter,:]*1e6), label="Simulation")
        sinusoidal_axes[j].plot(time_results, (sinusoid_data_harmonics[harm_counter,:]*1e6), label="Data")
        ax2=ramped_axes[j].twinx()
        ax3=sinusoidal_axes[j].twinx()
        ax2.set_yticks([])
        ax3.set_yticks([])
        ax2.set_ylabel(str(start_harm+harm_counter), rotation=0)
        ax3.set_ylabel(str(start_harm+harm_counter), rotation=0)
        if harm_counter == (len(ramped_harmonics)-1)/2:
            ramped_axes[j].set_ylabel("Current($\mu$A)")
            sinusoidal_axes[j].set_ylabel("Current($\mu$A)")
            ramped_axes[j].legend()
            sinusoidal_axes[j].legend()
        if harm_counter!=(len(ramped_harmonics)-1):
            ramped_axes[j].set_xticks([])
            sinusoidal_axes[j].set_xticks([])
        else:
            ramped_axes[j].set_xlabel("Time(s)")
            sinusoidal_axes[j].set_xlabel("Time(s)")

        j+=1
    #noramp_results.simulation_options["dispersion_bins"]=32


    #test_voltages=np.interp(noramp_results.t_nondim(noramp_results.time_vec[noramp_results.time_idx:]), time_results, voltage_results)
    #plt.plot(test_voltages,noramp_results.i_nondim(cmaes_time))
    #plt.plot(voltage_results[noramp_results.time_idx:], (current_results[noramp_results.time_idx:]), alpha=0.5)
    counter+=1
plt.show()
