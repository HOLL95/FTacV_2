import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_2"
file="Noramp_3_cv_fixed_ru.fixed_alpha"

master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","cap_phase","phase","alpha"]
noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
noramp_results.def_optim_list(master_optim_list)
method="timeseries"

cmaes_time=noramp_results.test_vals(param_vals, method)
noramp_results.simulation_options["likelihood"]=method
current_results=noramp_results.other_values["experiment_current"]
voltage_results=noramp_results.other_values["experiment_voltage"]
time_results=noramp_results.other_values["experiment_time"]
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
updated_lb=np.append(np.multiply(param_vals, 0.75),0.75*error)
updated_ub=np.append(np.multiply(param_vals, 1.25), 1.25*error)
for i in range(0, len(updated_lb)-1):
    print updated_lb[i],param_vals[i], updated_ub[i]
updated_b=[updated_lb, updated_ub]
updated_b=np.sort(updated_b, axis=0)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
print updated_b
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(param_vals, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
noramp_results.simulation_options["label"]="MCMC"
for i in range(0,8):
    mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
    mcmc.set_parallel(True)
    mcmc.set_max_iterations(50000)
    chains=mcmc.run()
    #pints.plot.trace(chains)
    #print file
    #plt.show()
    save_file="Noramp_3_"+method+"MCMC_"+str(i+1)+"_run2"
    f=open(("/").join([dir_path, results_dict,Electrode, "MCMC_runs", save_file]), "w")
    np.save(f, chains)
    f.close()
