import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
import math
import scipy.stats as stat
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_3"

final_osc=[25, 20, 15, 10]
other_files=["6", "9"]
for i in range(2, 11):

    if str(i) in other_files:
        file="Noramp_"+str(i)+"_cv_high_ru.run3_4"
    else:
        file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
    method="timeseries"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase", "phase", "alpha"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
    param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    #noramp_results.dim_dict["v_nondim"]=True
    #param_vals[2]=150
    noramp_results.def_optim_list(master_optim_list)
    #noramp_results.simulation_options["label"]="MCMC"
    noramp_results.def_optim_list(master_optim_list)
    print(param_vals)
    cmaes_time=noramp_results.test_vals(param_vals, method)

    #noramp_results.simulation_options["likelihood"]="fourier"
    noramp_results.simulation_options["dispersion_bins"]=16

    #reduced_idx=tuple(np.where(noramp_results.other_values["experiment_time"]<final_osc[i]))
    #cmaes_time=cmaes_time[reduced_idx]
    current_results=noramp_results.other_values["experiment_current"]#[reduced_idx]
    voltage_results=noramp_results.other_values["experiment_voltage"]#[reduced_idx]
    time_results=noramp_results.other_values["experiment_time"]#[reduced_idx]

    print(len(current_results))
    #noramp_results.time_vec=noramp_results.time_vec[np.where(noramp_results.time_vec<final_osc[i])]
    #print(len(current_results))



    if method=="timeseries":
        fit_data=current_results
        fit_times=time_results
    elif method=="fourier":
        fit_data=noramp_results.kaiser_filter(current_results)
        fit_times=np.linspace(0, 1, len(fit_data))
    error=np.std(np.abs(np.subtract(cmaes_time, fit_data)))
    mcmc_problem=pints.SingleOutputProblem(noramp_results, fit_times, fit_data)

    #updated_lb=np.append([noramp_results.param_bounds[key][0] for key in master_optim_list],0.75*error)
    #updated_ub=np.append([noramp_results.param_bounds[key][1] for key in master_optim_list], 1.25*error)
    error=np.std(abs(np.subtract(cmaes_time, current_results)))
    #plt.plot(voltage_results, cmaes_time)
    #plt.plot(voltage_results, current_results, alpha=0.5)
    #plt.show()
    mcmc_problem=pints.SingleOutputProblem(noramp_results, time_results, current_results)
    noramp_results.param_bounds["Ru"]=[0, 2000]
    noramp_results.param_bounds["k_0"]=[0, 250]
    noramp_results.param_bounds["alpha"]=[0.35, 0.615]
    updated_lb=np.append([noramp_results.param_bounds[key][0] for key in master_optim_list], 0.01*error)
    updated_ub=np.append([noramp_results.param_bounds[key][1] for key in master_optim_list], 10*error)
    [print(x, y, z) for x, y, z in zip(updated_lb, param_vals, updated_ub)]
    #updated_lb=np.append([x*0.65 for x in param_vals],0.01*error)
    #updated_ub=np.append([x*1.35 for x in param_vals], 10*error)
    updated_b=[updated_lb, updated_ub]
    updated_b=np.sort(updated_b, axis=0)
    #error=1
    log_liklihood=pints.GaussianLogLikelihood(mcmc_problem)
    #log_liklihood=pints.GaussianKnownSigmaLogLikelihood(mcmc_problem, error)
    #print(noramp_results.n_parameters(), len(updated_b[0]))
    log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
    #print(log_prior.n_parameters(), log_liklihood.n_parameters())
    log_posterior=pints.LogPosterior(log_liklihood, log_prior)
    #[(noramp_results.param_bounds[x][1]+noramp_results.param_bounds[x][0])/2 for x in noramp_results.optim_list ]
    mcmc_parameters=param_vals
    mcmc_parameters=np.append(mcmc_parameters, error)
    xs=[mcmc_parameters,
        mcmc_parameters,
        mcmc_parameters
        ]
    noramp_results.simulation_options["label"]="MCMC"
    noramp_results.simulation_options["test"]=False
    num_runs=20
    scores=np.ones(num_runs)*10
    skews=np.ones(num_runs)*10
    for j in range(0, num_runs):
        current_min=min(scores)
        mcmc = pints.MCMCController(log_posterior, 3, xs,method=pints.HaarioBardenetACMC)
        alpha_index=noramp_results.optim_list.index("alpha")
        alpha_chain=[]
        mcmc.set_parallel(False)
        mcmc.set_max_iterations(60000)
        chains=mcmc.run()
        rhat_mean=np.mean(pints.rhat_all_params(chains[:, 30000:, :]))
        for q in range(0, 2):
            alpha_chain=np.append(alpha_chain, chains[q, 30000:, alpha_index])
        alpha_skew=stat.skew(alpha_chain)
        if alpha_skew<-0.05:
            Electrode_save="Yellow"
            run2="MCMC_runs/omega_nondim/high_skew"
            save_file=file+"_MCMC_run9"
            filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run2])
            if abs(alpha_skew)<min([abs(x) for x in skews]) and rhat_mean<1.1:
                f=open(filepath+"/"+save_file, "wb")
                np.save(f, chains)
                f.close()
        else:
            Electrode_save="Yellow"
            run2="MCMC_runs/omega_nondim"
            save_file=file+"_MCMC_run9"
            filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run2])
        #print(pints.rhat_all_params(chains[:, 20000:, :]))
        #k_rhat=pints.rhat_all_params(chains[:, 20000:, :])[2]
        #pints.plot.trace(chains)
        #plt.show()
            if rhat_mean<1.08:
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
