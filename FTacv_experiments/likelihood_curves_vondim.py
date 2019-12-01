import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
import math
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_3"
def log_liklihood(sim,data, noise):
    n_t=len(data)
    distance=np.sum(np.power(np.subtract(data, sim),2))
    return -((n_t/2)*np.log(2*math.pi)) - (n_t*np.log(noise))-((1/2*noise**2)*distance)


file="Noramp_3_cv_high_ru.run3_2"
file2=os.getcwd()+"/Inferred_params/Yellow/MCMC_runs/omega_nondim"+"/Noramp_1_cv_high_ru_MCMC.run3"
chains=np.load(file2)
noise=np.mean(chains[:, 25000:, -1])
method="timeseries"
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])

noramp_results.def_optim_list(master_optim_list)
#noramp_results.simulation_options["label"]="MCMC"
noramp_results.def_optim_list(master_optim_list)
noramp_results.param_bounds["alpha"]=[0.3, 0.7]
noramp_results.param_bounds["Ru"]=[0, 1e4]

noramp_results.def_optim_list(noramp_results.optim_list)

normed_results=noramp_results.change_norm_group(param_vals, "norm")
cmaes_time=noramp_results.test_vals(param_vals, method)
current_results=noramp_results.other_values["experiment_current"]
voltage_results=noramp_results.other_values["experiment_voltage"]
time_results=noramp_results.other_values["experiment_time"]
noramp_results.simulation_options["no_transient"]=False
dim_time=noramp_results.t_nondim(time_results)
dim_voltage=noramp_results.e_nondim(voltage_results)
dim_current=noramp_results.i_nondim(current_results)
noramp_results.dim_dict["v_nondim"]=True
time_series=(noramp_results.test_vals(param_vals, method))
noramp_results.dim_dict["time_end"]=dim_time[-1]/noramp_results.nd_param.c_T0
#noramp_results.nd_param.time_end=time_results[-1]
#noramp_results.simulation_options["interpolant"]=time_results
#srs=[20, 50, 100, 200, 400, 1000, 2000]
#sfs=[1/x for x in srs]
#for j in range(0, len(sfs)):



noramp_results.times(1/5000.0)
print(noramp_results.nd_param.sampling_freq, len(noramp_results.time_vec))
time_series=(noramp_results.test_vals(param_vals, method))

#plt.plot(current_results)

#noramp_results.simulation_options["likelihood"]=method
#noramp_results.simulation_optio    ns["dispersion_bins"]=16
line_length=50
line_distance=0.05

increment1=np.linspace(-line_distance, 0, line_length+1)
increment2=np.linspace(0, line_distance, line_length-1)
increment=np.append(increment1, increment2[1:])
results=noramp_results.i_nondim(current_results)
scores=np.zeros(len(increment))
srs=[10, 25, 50, 100,200, 400, 800, 1000]#, 2000]
sfs=[1/x for x in srs]#, 200, 400, 800, 1000]]
#plt.plot(time_results, voltage_results)
#plt.show()
noramp_results.nd_param.time_end=time_results[-1]
for q in range(0, len(sfs)):
    #print(noramp_results.dim_dict["v_nondim"])
    noramp_results.times(sfs[q])
    for j in range(0, len(increment)):
        vals=noramp_results.change_norm_group(np.add(normed_results,increment[j]), "un_norm")
        time_series=(noramp_results.test_vals(vals, method))
        plot_current=np.interp(dim_time, noramp_results.t_nondim(noramp_results.time_vec), noramp_results.i_nondim(time_series))

        scores[j]=log_liklihood(dim_current,plot_current,noise)
    print(len(time_series))
    max_val=np.where(scores==max(scores))
    best_increase=increment[tuple(max_val)]
    if len(best_increase)>1:
        best_increase=best_increase[0]
    vals=noramp_results.change_norm_group(np.add(normed_results,best_increase), "un_norm")
    print(list(vals))
    time_series=noramp_results.i_nondim(noramp_results.test_vals(vals, method))
    #plt.plot(noramp_results.time_vec[noramp_results.time_idx:], time_series)
    #plt.plot(noramp_results.time_vec[noramp_results.time_idx:], interped_results)
    #plt.show()
    plt.subplot(2, 4, q+1)
    plt.plot(increment, scores, label=sfs[q])
    plt.axvline(best_increase, color="black", linestyle="--")
    plt.ylabel("Log-likelihood")
    plt.xlabel("Increment")
    plt.title("Sampling rate="+str(1/sfs[q]))
plt.show()
