import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
from params_class import params
import math
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_3"
def log_liklihood(sim,data, noise):
    n_t=len(data)
    distance=np.sum(np.power(np.subtract(data, sim),2))
    return -((n_t/2)*np.log(2*math.pi)) - (n_t*np.log(noise))-((1/2*noise**2)*distance)

for i in range(4, 7):
    file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
    file2="Noramp_"+str(i)+"_cv_high_ru_MCMC.run3"
    method="timeseries"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
    chain=np.load(("/").join([dir_path, results_dict,Electrode, "MCMC_runs","omega_nondim", file2]))
    noise=np.mean(chain[:, 25000:, -1])
    param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    noramp_results.def_optim_list(master_optim_list)
    #noramp_results.simulation_options["label"]="MCMC"
    noramp_results.def_optim_list(master_optim_list)
    noramp_results.param_bounds["alpha"]=[0.3, 0.7]
    noramp_results.param_bounds["Ru"]=[0, 1e4]

    noramp_results.def_optim_list(noramp_results.optim_list)

    param_vals_2=[0.233, 0.073, 343.621, 65, 7.602e-5, 8.026e-3, -6.688e-4, 6.445e-11, 8.941, 4.211,4.68, 0.65]
    normed_results=noramp_results.change_norm_group(param_vals, "norm")
    normed_results_2=noramp_results.change_norm_group(param_vals_2, "norm")
    line_equation=[x-y for x,y in zip(normed_results, normed_results_2)]



    cmaes_time=noramp_results.test_vals(param_vals, method)
    current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])
    voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])
    time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])
    noramp_results.nd_param.time_end=time_results[-1]
    noramp_results.simulation_options["interpolant"]=time_results
    #srs=[20, 50, 100, 200, 400, 1000, 2000]
    #sfs=[1/x for x in srs]
    #for j in range(0, len(sfs)):

    time_series=noramp_results.i_nondim(noramp_results.test_vals(param_vals_2, method))
    line_length=50
    end_length=5
    line_distance=0.05
    total_line=np.linspace(-0.1,1.1 , line_length)
    increment=np.diff(total_line)[0]
    noramp_results.dim_dict["v_nondim"]=True
    noramp_results.simulation_options["no_transient"]=False
    noramp_results.nd_param=params(noramp_results.dim_dict)
    #total_line=np.append(total_line, positive_end)
    #total_line=line
    param_list=[noramp_results.change_norm_group(np.add(normed_results_2, np.multiply(x, line_equation)), "un_norm") for x in total_line]
    #param_list=[np.add(normed_results_2, np.multiply(x, line_equation)) for x in total_line]
    #noramp_results.simulation_options["likelihood"]=method
    #noramp_results.simulation_optio    ns["dispersion_bins"]=16

    scores=np.zeros(len(total_line))
    sfs=[1/x for x in [20, 50, 100, 200, 500, 1000]]#, 200, 400, 800, 1000]]
    #plt.plot(time_results, voltage_results)
    #plt.show()

    noramp_results.nd_param.time_end=time_results[-1]
    for q in range(0, len(sfs)):
        #print(noramp_results.dim_dict["v_nondim"])
        noramp_results.dim_dict["sampling_freq"]=sfs[q]
        noramp_results.dim_dict["time_end"]=time_results[-1]
        noramp_results.nd_param=params(noramp_results.dim_dict)
        noramp_results.times()
        print(noramp_results.nd_param.c_T0)
        for j in range(0, len(total_line)):
            #vals=noramp_results.change_norm_group(np.add(normed_results,increment[j]), "un_norm")
            time_series=noramp_results.i_nondim(noramp_results.test_vals(param_list[j], method))
            syn_times=noramp_results.t_nondim(noramp_results.time_vec)
            interped_results=np.interp(time_results,syn_times,time_series)
            #plt.plot(voltage_results,interped_results)
            #plt.plot(voltage_results, current_results)
            #plt.plot(time_results, current_results)
            #plt.show()
            #plt.show()
            scores[j]=log_liklihood(current_results,interped_results,noise)
        print(len(time_series))
        #max_val=np.where(scores==max(scores))
        #best_increase=increment[tuple(max_val)]
        #if len(best_increase)>1:
        #    best_increase=best_increase[0]
        #vals=noramp_results.change_norm_group(np.add(normed_results,best_increase), "un_norm")
        #time_series=noramp_results.i_nondim(noramp_results.test_vals(vals, method))
        #plt.plot(noramp_results.time_vec[noramp_results.time_idx:], time_series)
        #plt.plot(noramp_results.time_vec[noramp_results.time_idx:], interped_results)
        #plt.show()
        plt.subplot(2, 3, q+1)
        plt.plot(total_line, [noramp_results.normalise(x, [min(scores), max(scores)]) for x in scores], label=sfs[q])
        #plt.axvline(best_increase, color="black", linestyle="--")
        plt.ylabel("Relative log-likelihood")
        plt.xlabel("Increment")
        plt.title("Sampling rate="+str(1/sfs[q]))
    plt.show()
