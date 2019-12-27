import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
from harmonics_plotter import harmonics
import pints.plot
import os
from multiplotter import multiplot
import math
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
ramped_data_class.simulation_options["dispersion_bins"]=50
ramped_data_path=("/").join([dir_path, "experiment_data_2", "Experimental-120919", "Ramped", "Yellow"])
start_harm=2
end_harm=7
ramp_data_harm_class=harmonics(list(range(start_harm, end_harm)), ramped_data_class.nd_param.omega, 0.05)
harmonic_range=list(range(start_harm, end_harm))
ramped_data_class.harmonic_range=harmonic_range
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
ramped_current_results=np.loadtxt(ramped_data_path+"/Yellow_Electrode_Ramped_"+str(i)+"_cv_current")
ramped_voltage_results=np.loadtxt(ramped_data_path+"/Yellow_Electrode_Ramped_"+str(i)+"_cv_voltage")[:,1]

ramped_time_results=ramped_current_results[:,0]
ramped_current_results=ramped_current_results[:,1]
ramped_data_class.other_values["experiment_time"]=ramped_time_results
ramped_data_class.other_values["experiment_current"]=ramped_current_results
ramped_data_class.other_values["experiment_voltage"]=ramped_voltage_results
ramp_fit=single_electron(None, dim_paramater_dictionary=ramped_data_class.dim_dict, simulation_options=ramped_data_class.simulation_options, other_values=ramped_data_class.other_values, param_bounds=ramped_data_class.param_bounds)
ramped_data_harmonics=ramp_data_harm_class.generate_harmonics(ramped_time_results, ramped_current_results)
file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","phase", "alpha"]
ramp_fit.def_optim_list(master_optim_list)
ramp_fit=single_electron(("/").join([dir_path, results_dict,Electrode,"Run_3", file]))
noramp_harm=harmonics(list(range(start_harm, end_harm)), ramp_fit.nd_param.omega, 0.05)
param_vals=([ramp_fit.save_dict["params"][0][ramp_fit.save_dict["optim_list"].index(key)] if  (key in ramp_fit.save_dict["optim_list"]) else ramp_fit.dim_dict[key] for key in master_optim_list])
ramp_fit.def_optim_list(master_optim_list)

#print(param_vals)
cmaes_time_series=ramp_fit.i_nondim(ramp_fit.test_vals(param_vals, method))
current_results=ramp_fit.i_nondim(ramp_fit.other_values["experiment_current"])
voltage_results=ramp_fit.e_nondim(ramp_fit.other_values["experiment_voltage"])
time_results=ramp_fit.t_nondim(ramp_fit.other_values["experiment_time"])
cmaes_times=time_results#ramp_fit.t_nondim(ramp_fit.time_vec[ramp_fit.time_idx:])
test_voltages=voltage_results
filtered_exp, filtered_exp_freq=likelihood_func(harmonic_range, cmaes_times, cmaes_time_series, ramp_fit.nd_param.omega)
filtered_data, filtered_data_freq=likelihood_func(harmonic_range, time_results, current_results,  ramp_fit.nd_param.omega)

#test_voltages=np.interp(cmaes_times, time_results, voltage_results)
ramped_param_vals=[ramp_fit.dim_dict[x] for x in ramped_optim_list]
ramped_param_vals=[0.230288480903605996, 0.041226544633855209, 122.16602260400302, 742.0285824092344, 8.212777750912338e-05, 0.00295011584929738, -0.0005139373308230138, 7.330227126430568e-11, 8.885659818458482,0,0.5146469208293426]
values=[[0.23794865088573297, 0.012233237509352422, 127.31889407388738, 866.5594838874783, 7.735179022904297e-05, 0.002100658721047255, -0.0003313772993110888, 7.92266323442856e-11, 8.884799587792013, 0, 0.5999996893552537],
        [0.23363378517047495, 0.03481010462713212, 125.2418680556816, 866.5594838874783, 7.694171200331618e-05, 0.003209611999861764, -0.0004263185805571494, 7.476933579891946e-11, 8.884799587792013, 5.02832295276801, 0.5999996422725197] ,
        [0.23422376456952138, 0.03260990314447788, 127.84406477439944, 866.5594838874783, 7.6536077935438e-05, 0.0028147360128534457, -0.00040069155469661145, 7.510752378483546e-11, 8.884799587792013,  5.040933894042927, 0.5999999355676289]]
#axes=multiplot(1, len(values)-1, **{"harmonic_position":0, "num_harmonics":5, "orientation":"landscape", "plot_width":5})

#ramped_param_vals=[0.2146905833892083, 0.03054987719700828, 97.73281867896537, 607.3604383735801, 9.855333151006176e-05, 0.002360092744665423, -0.0006167247956683124, 2.816480836302093e-10, 8.884799587792013, 3.7667007049769547, 5.036133769048488, 0.6175763030220238]
#ramped_param_vals=[0.2105491969962405, 0.060492147323936776, 109.9494203694055, 667.8257252271911, 9.034055522879727e-05, 0.002655104265433323, -0.0005653310639016844, 1.1167514700671097e-10, 8.884680131969208, 1.232169931719665, 0.46318222889196947]
#ramped_param_vals=[0.228308296004341, 0.04678928704199179, 109.9494215030595, 753.3526020630202, 7.39149998133382e-05, 0.0032451274206835596, -0.00046254359814023406, 8.063249835588253e-11, 8.884745411351478, 2.236534676329661, 0.4631822303721044]
#ramped_param_vals=[0.22285754799270494, 0.04694179204711807, 247.23685999878182, 512.1240571770394, 7.706083559507504e-05, 0.0028478498239308467, -0.00041016659783854256, 7.257026464642658e-11, 8.885659818458482,0, 0.6999999975394792]
"""
ramp_fit.harmonic_range=harmonic_range
axes_keys=sorted(axes.axes_dict.keys())
j=0

for i in range(1, len(values)):
    ramped_time_series=ramp_fit.i_nondim(ramp_fit.test_vals(values[i], "timeseries"))
    ramped_times=ramp_fit.t_nondim(ramp_fit.time_vec[ramp_fit.time_idx:])
    ramped_harmonics=ramp_data_harm_class.generate_harmonics(ramped_times, ramped_time_series)
    results=np.loadtxt(dir_path+"/experiment_data_2/Experimental-120919/Ramped/Yellow/"+"Yellow_Electrode_Ramped_"+str(i+1)+"_cv_current")
    time_results=results[:,0]
    current_results=results[:,1]
    sinusoid_data_harmonics=noramp_harm.generate_harmonics(time_results, current_results)
    sinusoid_harmonics=noramp_harm.generate_harmonics(cmaes_times, cmaes_time_series)
    ramp_fit.simulation_options["method"]="dcv"
    dcv_plot=ramp_fit.e_nondim(ramp_fit.define_voltages())
    ramp_fit.simulation_options["method"]="ramped"
    dcv_data_plot=np.interp(ramped_time_results, ramped_times, dcv_plot)
    for harm_counter in range(0, len(ramped_harmonics)):
        ax=axes.axes_dict["row1"][j]
        ax.plot(dcv_plot, (ramped_harmonics[harm_counter,:]*1e6), label="Simultation")
        ax.plot(dcv_data_plot, (ramped_data_harmonics[harm_counter,:]*1e6), label="Data", alpha=0.5)
        ax.axvline(values[i][0], linestyle="--", color="black")
        ax2=ax.twinx()
        ax2.set_ylabel(harmonic_range[harm_counter], rotation=0)
        ax2.set_yticks([])
            #plt.axvline(ramped_param_vals[0], color="black", linestyle="--", label="$E_0$")
        if harm_counter%3==2:
            ax.set_ylabel("Current($\mu$A)")
        if harm_counter==(len(ramped_harmonics)-1):
            ax.set_xlabel("Time(s)")
            ax.legend()
        else:
            ax.set_xticks([])
        #sinusoidal_axes[j].plot(cmaes_times, (sinusoid_harmonics[harm_counter,:]*1e6), label="Simulation")
        #sinusoidal_axes[j].plot(time_results, (sinusoid_data_harmonics[harm_counter,:]*1e6), label="Data")
        j+=1



plt.show()
"""
ramped_data_class.harmonic_range=range(3, 7)
print(len(ramped_data_class.time_vec), len(ramped_time_results))
ramped_data_class.def_optim_list(master_optim_list)
fourier_arg=ramped_data_class.kaiser_filter(ramped_data_class.other_values["experiment_current"])
ramped_data_class.time_vec=ramped_data_class.other_values["experiment_time"]
ramped_time_series=ramped_data_class.test_vals(values[-1], "timeseries")
fourier_test=ramped_data_class.kaiser_filter(ramped_time_series)
error=np.std(abs(np.subtract(fourier_test, fourier_arg)))

ramped_data_class.simulation_options["method"]="ramped"
ramped_data_class.simulation_options["likelihood"]="fourier"
"""
ramped_data_class.def_optim_list(["E0_mean","k_0","Ru","Cdl","CdlE1", "CdlE2","omega","gamma", "alpha"])
for i in range(0, len(ramped_data_class.optim_list)):
    idx=master_optim_list.index(ramped_data_class.optim_list[i])
    ramped_data_class.param_bounds[ramped_data_class.optim_list[i]]=[0.9*ramped_param_vals[idx], 1.1*ramped_param_vals[idx]]
"""
print(ramped_data_class.param_bounds)
ramped_data_class.param_bounds["Ru"]=[0, 1e3]
ramped_data_class.param_bounds["k_0"]=[0, 250]
ramped_data_class.param_bounds["alpha"]=[0.4, 0.6]
ramped_data_class.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","phase", "alpha"])
ramped_data_class.dim_dict["cap_phase"]=0
#ramped_param_vals=[0.22788480903605996, 0.1326544633855209, 122.16602260400302, 742.0285824092344, 8.212777750912338e-05, 0.00295011584929738, -0.0005139373308230138, 7.130227126430568e-11, 8.885659818458482,0,0,0.5146469208293426]
mcmc_problem=pints.SingleOutputProblem(ramped_data_class, np.linspace(0, 1, len(fourier_test)), fourier_arg)
ramped_data_class.secret_data_fourier=fourier_arg
ramped_data_class.simulation_options["label"]="MCMC"
#ramped_data_class.simulation_options["test"]=True
num_runs=10
updated_lb=np.append([ramped_data_class.param_bounds[key][0] for key in master_optim_list], 0.01*error)
updated_ub=np.append([ramped_data_class.param_bounds[key][1] for key in master_optim_list], 10*error)
[print(x, y, z) for x, y, z in zip(updated_lb, param_vals, updated_ub)]
#updated_lb=np.append([x*0.65 for x in param_vals],0.01*error)
#updated_ub=np.append([x*1.35 for x in param_vals], 10*error)
updated_b=[updated_lb, updated_ub]
updated_b=np.sort(updated_b, axis=0)
log_liklihood=pints.GaussianLogLikelihood(mcmc_problem)
#log_liklihood=pints.GaussianKnownSigmaLogLikelihood(mcmc_problem, error)
#print(ramped_data_class.n_parameters(), len(updated_b[0]))
log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
#print(log_prior.n_parameters(), log_liklihood.n_parameters())
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
#[(ramped_data_class.param_bounds[x][1]+ramped_data_class.param_bounds[x][0])/2 for x in ramped_data_class.optim_list ]
mcmc_parameters=values[-1]
mcmc_parameters=np.append(mcmc_parameters, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
ramped_data_class.simulation_options["label"]="MCMC"
ramped_data_class.simulation_options["test"]=False
num_runs=20
scores=np.ones(num_runs)*10
skews=np.ones(num_runs)*10
for j in range(0, num_runs):
    current_min=min(scores)
    mcmc = pints.MCMCController(log_posterior, 3, xs,method=pints.HaarioBardenetACMC)
    alpha_index=ramped_data_class.optim_list.index("alpha")
    alpha_chain=[]
    mcmc.set_parallel(False)
    mcmc.set_max_iterations(100000)
    chains=mcmc.run()
    rhat_mean=np.mean(pints.rhat_all_params(chains[:, 60000:, :]))
    for q in range(0, 2):
        alpha_chain=np.append(alpha_chain, chains[q, 60000:, alpha_index])
    alpha_skew=stat.skew(alpha_chain)
    if alpha_skew<-0.05:
        Electrode_save="Yellow/Ramped"
        run2="MCMC_runs/v_nondim/high_skew"
        save_file=ramped_file+"_MCMC_run8"
        filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run2])
        if abs(alpha_skew)<min([abs(x) for x in skews]):
            f=open(filepath+"/"+save_file, "wb")
            np.save(f, chains)
            f.close()
    else:
        Electrode_save="Yellow/Ramped"
        run2="MCMC_runs/v_nondim"
        save_file=ramped_file+"_MCMC_run8"
        filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run2])
    #print(pints.rhat_all_params(chains[:, 20000:, :]))
    #k_rhat=pints.rhat_all_params(chains[:, 20000:, :])[2]
    #pints.plot.trace(chains)
    #plt.show()
        if rhat_mean<1.1:
            f=open(filepath+"/"+save_file, "wb")
            np.save(f, chains)
            f.close()
            break
        elif rhat_mean<min(scores):
            f=open(filepath+"/"+save_file, "wb")
            np.save(f, chains)
            f.close()
    scores[j]=rhat_mean
    skews[j]=alpha_skew
