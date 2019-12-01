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
ramp_data_class=single_electron(("/").join([dir_path, results_dict,Electrode,run, ramped_file]), {}, {}, {}, {}, False)
ramped_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
ramp_data_class.def_optim_list(ramped_optim_list)
ramp_data_class.dim_dict["v_nondim"]=True
ramp_data_class.dim_dict["CdlE3"]=0
ramp_data_class.dim_dict["phase"]=0
ramp_data_class.dim_dict["cap_phase"]=0
#ramp_data_class.dim_dict["omega"]=8.95
ramp_data_class.dim_dict["original_omega"]=8.94
ramp_data_class.simulation_options["dispersion_bins"]=16
ramped_data_path=("/").join([dir_path, "experiment_data_2", "Experimental-221119", "Ramped", "Red"])
start_harm=2
end_harm=6
ramp_data_harm_class=harmonics(list(range(start_harm, end_harm)), ramp_data_class.nd_param.omega, 0.1)
harmonic_range=list(range(start_harm, end_harm))
ramp_data_class.harmonic_range=harmonic_range
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

i="2"
ramped_current_results=np.loadtxt(ramped_data_path+"/Ramped_red_1e0M_scan1_cv_current")
ramped_voltage_results=np.loadtxt(ramped_data_path+"/Ramped_red_1e0M_scan1_cv_voltage")[:,1]

ramped_time_results=ramped_current_results[:,0]
ramped_current_results=ramped_current_results[:,1]
ramp_data_class.other_values["experiment_time"]=ramped_time_results
ramp_data_class.other_values["experiment_current"]=ramped_current_results
ramp_data_class.other_values["experiment_voltage"]=ramped_voltage_results
ramp_fit=single_electron(None, dim_paramater_dictionary=ramp_data_class.dim_dict, simulation_options=ramp_data_class.simulation_options, other_values=ramp_data_class.other_values, param_bounds=ramp_data_class.param_bounds)
ramped_data_harmonics=ramp_data_harm_class.generate_harmonics(ramped_time_results, ramped_current_results)
file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","phase", "alpha"]
ramp_fit.def_optim_list(master_optim_list)
noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,"Run_3", file]))
noramp_harm=harmonics(list(range(start_harm, end_harm)), noramp_results.nd_param.omega, 0.5)
param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
noramp_results.def_optim_list(master_optim_list)

#print(param_vals)
cmaes_time_series=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))
current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])
voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])
time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])
cmaes_times=time_results#noramp_results.t_nondim(noramp_results.time_vec[noramp_results.time_idx:])
test_voltages=voltage_results


#test_voltages=np.interp(cmaes_times, time_results, voltage_results)
ramped_param_vals=[noramp_results.dim_dict[x] for x in ramped_optim_list]
#ramped_param_vals=[0.19148814061983946, 0.03702918672685089, 126.61824647795085, 1.0128942825395, 8.296772457318047e-05, 0.010335675369948137, 0.0032277396415101856, 1.60187250145465e-10, 8.942010060550901, 0, 0.5]
#ramped_param_vals=[0.21803890811600762, 0.14559210666411324, 59.49933750290762, 845.8134700849996, 3.856873683824087e-05, 0.020501659307914116, -8.273937707160492e-05, 6.416228744721123e-11, 8.94240058250762, 0, 0.30000000260617204]
ramped_param_vals=[0.16569534640804542, 0.07406404751300973, 74.45004247951917, 100.2855105517765, 6.982979823375169e-05, 0.07580732986337356, 0.0016387886539689818, 6.498771582276014e-10, 8.941802744833884, 0, 0.300001923250468]


#ramped_param_vals=[0.230288480903605996, 0.041226544633855209, 122.16602260400302, 742.0285824092344, 8.212777750912338e-05, 0.00295011584929738, -0.0005139373308230138, 7.330227126430568e-11, 8.945659818458482,0,0.5146469208293426]
#ramped_param_vals=[0.2146905833892083, 0.03054987719700828, 97.73281867896537, 607.3604383735801, 9.855333151006176e-05, 0.002360092744665423, -0.0006167247956683124, 2.816480836302093e-10, 8.884799587792013, 3.7667007049769547, 5.036133769048488, 0.6175763030220238]
#ramped_param_vals=[0.2105491969962405, 0.060492147323936776, 109.9494203694055, 667.8257252271911, 9.034055522879727e-05, 0.002655104265433323, -0.0005653310639016844, 1.1167514700671097e-10, 8.884680131969208, 1.232169931719665, 0.46318222889196947]
#ramped_param_vals=[0.228308296004341, 0.04678928704199179, 109.9494215030595, 753.3526020630202, 7.39149998133382e-05, 0.0032451274206835596, -0.00046254359814023406, 8.063249835588253e-11, 8.884745411351478, 2.236534676329661, 0.4631822303721044]

ramp_fit.harmonic_range=harmonic_range
e_s=[ramped_param_vals[0]]#[0.222, 0.228, 0.234, 0.236]
for i in e_s:
    ramped_param_vals[0]=i
    ramped_time_series=ramp_fit.i_nondim(ramp_fit.test_vals(ramped_param_vals, "timeseries"))
    #plt.plot(ramped_time_series)
    #plt.show()
    ramped_times=ramp_fit.t_nondim(ramp_fit.time_vec[ramp_fit.time_idx:])
    ramped_harmonics=ramp_data_harm_class.generate_harmonics(ramped_times, ramped_time_series)
    sinusoid_data_harmonics=noramp_harm.generate_harmonics(time_results, current_results)
    sinusoid_harmonics=noramp_harm.generate_harmonics(cmaes_times, cmaes_time_series)
    ramp_fit.simulation_options["method"]="dcv"
    dcv_plot=ramp_fit.e_nondim(ramp_fit.define_voltages())
    ramp_fit.simulation_options["method"]="ramped"
    j=0
    for harm_counter in range(0, len(ramped_harmonics)):
        plt.subplot(len(ramped_harmonics), 1, j+1)
        plt.plot(dcv_plot, abs(ramped_harmonics[harm_counter,:]*1e6), label="Simultation")
        #if i==e_s[0]:
        plt.plot(dcv_plot, abs(ramped_data_harmonics[harm_counter,:]*1e6), label="Data")#, linestyle="--")
            #plt.axvline(ramped_param_vals[0], color="black", linestyle="--", label="$E_0$")
        plt.ylabel("Current($\mu$A)")
        #sinusoidal_axes[j].plot(cmaes_times, (sinusoid_harmonics[harm_counter,:]*1e6), label="Simulation")
        #sinusoidal_axes[j].plot(time_results, (sinusoid_data_harmonics[harm_counter,:]*1e6), label="Data")
        j+=1

    plt.xlabel("Time(s)")
    plt.legend()
plt.show()
