import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
file="Noramp_3_cv_high_ru.fourier"
save_file="Noramp_3_fourier_MCMC"
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","cap_phase","phase","alpha"]

noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode, file]))
param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
noramp_results.def_optim_list(master_optim_list)
cmaes_time=noramp_results.test_vals(param_vals, "timeseries")
noramp_results.simulation_options["likelihood"]="timeseries"
current_results=noramp_results.other_values["experiment_current"]
voltage_results=noramp_results.other_values["experiment_voltage"]
time_results=noramp_results.other_values["experiment_time"]
error=np.std(np.subtract(cmaes_time, current_results))
mcmc_problem=pints.SingleOutputProblem(noramp_results, time_results, current_results)
noramp_results.param_bounds["Ru"][1]=700
updated_lb=np.append([noramp_results.param_bounds[key][0] for key in master_optim_list],0)
updated_ub=np.append([noramp_results.param_bounds[key][1] for key in master_optim_list], 2*error)
for i in range(0, len(updated_lb)-1):
    print updated_lb[i],param_vals[i], updated_ub[i]
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_lb,
                                updated_ub)
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(param_vals, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
noramp_results.simulation_options["label"]="MCMC"
for i in range(0,5):
    mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
    mcmc.set_parallel(True)
    mcmc.set_max_iterations(10000)
    chains=mcmc.run()
pints.plot.trace(chains)
print file
plt.show()
f=open(("/").join([dir_path, results_dict,Electrode, save_file]), "w")
np.save(f, chains)
f.close()
