import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
import itertools
import math
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Red"
run="Run_1"
concs=["1e-1M", "1e0M"]
file_numbers=[str(x) for x in range(1, 4)]
list1=["Noramp__"+x+"_"+y+"_cv.run1" for x in concs for y in file_numbers]
file_list=['Noramp__1e-1M_1_cv.run1', 'Noramp__1e-1M_2_cv.run1', 'Noramp__1e0M_1_cv.run1', 'Noramp__1e0M_3_cv.run1']

print(file_list)
all_param_vals=[
    [0.19148814061983946, 0.03702918672685089, 126.61824647795085, 4402.0128942825395, 8.296772457318047e-05, 0.10335675369948137, 0.0032277396415101856, 9.60187250145465e-10, 8.942010060550901, 4.43785397688034, 6.143447241208483, 0.3],
    [0.18654805784839637, 0.03958691204050591, 108.44034734915691, 4509.974465151381, 7.597911463180463e-05, 0.10329644861871086, 0.0033984772908889297, 8.975717009117107e-10, 8.942056374476506, 4.42926027420097, 6.156393021351812, 0.3],
    [0.1857171783065358, 0.06967714147489301, 153.25311776193467, 3125.6277376474554, 9.999986745027279e-05, 0.054696607481221, -0.000210671324025637, 6.157332307641423e-10, 8.941736880130359, 4.4441854786799135, 1.6080016213805508e-07, 0.3],
    [0.1797867762426466, 0.07071588064415471, 121.54476416303484, 3284.1933603344582, 9.999997720150796e-05, 0.06008173605580561, 0.00031272007095725814, 7.228802704659992e-10, 8.94142857848522, 4.44186097843843, 7.182074176351618e-08, 0.3]
]
def likelihood_func(harmonic_range, times, time_series, omega):
    filter_val=0.5
    frequency=np.fft.fftfreq(len(times), times[1]-times[0])
    hanning=np.hanning(len(times))
    time_series=np.multiply(time_series, hanning)
    desired_idx=np.where((frequency>((harmonic_range[0]*omega)-(omega*filter_val)))&(frequency<(harmonic_range[-1]*omega)+(omega*filter_val)))
    Y=(np.fft.fft(time_series))
    return Y[desired_idx],frequency[desired_idx]
method="timeseries"
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
fig, ax=plt.subplots(len(file_list), 3)
for i in range(0, len(file_list)):
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file_list[i]]))
    #For 1e-1M no. 1
    #[0.19148814061983946, 0.03702918672685089, 126.61824647795085, 4402.0128942825395, 8.296772457318047e-05, 0.10335675369948137, 0.0032277396415101856, 9.60187250145465e-10, 8.942010060550901, 4.43785397688034, 6.143447241208483, 0.5]
    #For 1e-1M no. 2
    #[0.18654805784839637, 0.03958691204050591, 108.44034734915691, 4509.974465151381, 7.597911463180463e-05, 0.10329644861871086, 0.0033984772908889297, 8.975717009117107e-10, 8.942056374476506, 4.42926027420097, 6.156393021351812, 0.5]
    #For 1e0M no. 1
    #[0.1857171783065358, 0.06967714147489301, 153.25311776193467, 3125.6277376474554, 9.999986745027279e-05, 0.054696607481221, -0.000210671324025637, 6.157332307641423e-10, 8.941736880130359, 4.4441854786799135, 1.6080016213805508e-07, 0.30000003949418524]
    #For 1e0M no.3
    #[0.1797867762426466, 0.07071588064415471, 121.54476416303484, 3284.1933603344582, 9.999997720150796e-05, 0.06008173605580561, 0.00031272007095725814, 7.228802704659992e-10, 8.94142857848522, 4.44186097843843, 7.182074176351618e-08, 0.5]

    #param_vals=([noramp_results.save_dict["params"][1][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    param_vals=all_param_vals[i]
    noramp_results.def_optim_list(master_optim_list)
    #noramp_results.simulation_options["label"]="MCMC"
    cmaes_time=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))
    #noramp_results.simulation_options["likelihood"]=method
    #noramp_results.simulation_options["dispersion_bins"]=16
    dec_amount=16
    current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])
    voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])
    time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])
    harmonic_range=range(2, 8)
    fourier_sim, fourier_sim_freq=likelihood_func(harmonic_range, time_results, cmaes_time, noramp_results.dim_dict["omega"])
    data_sim, data_sim_freq=likelihood_func(harmonic_range, time_results, current_results, noramp_results.dim_dict["omega"])
    
    ax[i, 0].set_title(file_list[i])
    ax[i, 0].plot(voltage_results, (cmaes_time), label="Sim")
    ax[i,0].plot(voltage_results, ((current_results)), alpha=0.7, label="Data")
    ax[i, 0].set_xlabel("Voltage(V)")
    ax[i, 0].set_ylabel("Current(A)")
    ax[i, 0].legend()
    ax[i, 1].plot(fourier_sim_freq, np.real(fourier_sim), label="Sim")
    ax[i,1].plot(data_sim_freq, np.real(data_sim), alpha=0.7, label="Data")
    ax[i, 1].set_xlabel("Frequency(Hz)")
    ax[i, 1].set_ylabel("Real magnitude")
    ax[i, 1].legend()
    ax[i, 2].plot(fourier_sim_freq, np.imag(fourier_sim), label="Sim")
    ax[i, 2].set_xlabel("Frequency(Hz)")
    ax[i, 2].set_ylabel("Imaginary magnitude")
    ax[i,2].plot(data_sim_freq, np.imag(data_sim), alpha=0.7, label="Data")
    ax[i, 2].legend()
plt.show()
