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



for i in range(4, 7):
    file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
    file2="Noramp_"+str(i)+"_cv_high_ru_MCMC.run3"
    method="timeseries"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
    chain=np.load(("/").join([dir_path, results_dict,Electrode, "MCMC_runs","omega_nondim", file2]))
    noise=np.mean(chain[:, 25000:, -1])
    param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    #noramp_results.dim_dict["v_nondim"]=True
    noramp_results.def_optim_list(master_optim_list)
    #noramp_results.simulation_options["label"]="MCMC"
    noramp_results.def_optim_list(master_optim_list)
    print(param_vals)
    cmaes_time=noramp_results.test_vals(param_vals, method)
    #noramp_results.simulation_options["likelihood"]=method
    #noramp_results.simulation_options["dispersion_bins"]=16

    current_results=noramp_results.other_values["experiment_current"]
    voltage_results=noramp_results.other_values["experiment_voltage"]
    time_results=noramp_results.other_values["experiment_time"]

    plt.plot(noramp_results.t_nondim(time_results), current_results)
    plt.show()
    test_voltages=np.interp(noramp_results.t_nondim(noramp_results.time_vec[noramp_results.time_idx:]), time_results, voltage_results)
    plt.plot(noramp_results.t_nondim(noramp_results.time_vec[noramp_results.time_idx:]), test_voltages)
    plt.plot(time_results, voltage_results)
    plt.show()
    plt.plot(voltage_results,noramp_results.i_nondim(cmaes_time))
    #plt.plot(voltage_results[noramp_results.time_idx:], noramp_results.i_nondim((current_results[noramp_results.time_idx:])), alpha=0.5)
    plt.show()
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
mcmc_problem=pints.SingleOutputProblem(noramp_results, time_results, current_results)
updated_lb=[noramp_results.param_bounds[x][0] for x in noramp_results.optim_list]#np.append(,0.01*error)
updated_ub=[noramp_results.param_bounds[x][1] for x in noramp_results.optim_list]#np.append(, 10*error)
updated_b=[updated_lb, updated_ub]
updated_b=np.sort(updated_b, axis=0)
#error=1
log_liklihood=pints.GaussianKnownSigmaLogLikelihood(mcmc_problem, error)
#print(noramp_results.n_parameters(), len(updated_b[0]))
log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
#print(log_prior.n_parameters(), log_liklihood.n_parameters())
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
#[(noramp_results.param_bounds[x][1]+noramp_results.param_bounds[x][0])/2 for x in noramp_results.optim_list ]
mcmc_parameters=param_vals
#mcmc_parameters=np.append(mcmc_parameters, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
noramp_results.simulation_options["label"]="MCMC"
for i in range(0, 10):
    mcmc = pints.MCMCController(log_posterior, 3, xs,method=pints.HaarioBardenetACMC)
    mcmc.set_parallel(True)
    mcmc.set_max_iterations(20000)
    chains=mcmc.run()
    rhat_mean=np.mean(pints.rhat_all_params(chains))
    #pints.plot.trace(chains)
    #plt.show()
    if rhat_mean<1.2:
        Electrode="Yellow"
        Electrode_save="Yellow"
        folder="Noramp"
        Method="3_cv"
        sim_options="fixed_error_"+str(error)
        run="MCMC_runs/omega_nondim"
        save_file=filename=("_").join([folder,Method, sim_options, "MCMC"])+".run3_2"
        filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run])
        f=open(filepath+"/"+filename, "wb")
        np.save(f, chains)
        f.close()
        break
