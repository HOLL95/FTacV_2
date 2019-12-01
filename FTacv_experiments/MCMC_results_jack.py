import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
import math
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Red"
run="Run_1"


file="Noramp__1e0M_3_cv.run1"
method="timeseries"
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
#For 1e-1M no. 3
#noramp_results.save_dict["params"][2]=[0.1865595627995186, 0.03972665649139215, 108.05068250587402, 4508.805691705821, 7.59485472065489e-05, 0.10346815397552038, 0.0033306014497631285, 8.937287620154294e-10, 8.942057998043746, 4.429100575412366, 6.156899544154524, 0.30000014728427155]
param_vals=([noramp_results.save_dict["params"][2][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])

noramp_results.def_optim_list(master_optim_list)
#noramp_results.simulation_options["label"]="MCMC"
param_vals[-1]=0.5
print(param_vals)
cmaes_time=noramp_results.test_vals(param_vals, method)
#noramp_results.simulation_options["likelihood"]=method
#noramp_results.simulation_options["dispersion_bins"]=16
dec_amount=16
current_results=noramp_results.other_values["experiment_current"][0::dec_amount]
voltage_results=noramp_results.other_values["experiment_voltage"][0::dec_amount]
time_results=noramp_results.other_values["experiment_time"][0::dec_amount]
test_voltages=np.interp(noramp_results.t_nondim(noramp_results.time_vec[noramp_results.time_idx:]), time_results, voltage_results)
plt.plot(noramp_results.other_values["experiment_voltage"],noramp_results.i_nondim(cmaes_time))
plt.plot(voltage_results, noramp_results.i_nondim((current_results)), alpha=0.5)
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
updated_lb=np.append([noramp_results.param_bounds[x][0] for x in noramp_results.optim_list],0.01*error)
updated_ub=np.append([noramp_results.param_bounds[x][1] for x in noramp_results.optim_list], 10*error)
updated_b=[updated_lb, updated_ub]
updated_b=np.sort(updated_b, axis=0)
#error=1
log_liklihood=pints.GaussianLogLikelihood(mcmc_problem)
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
for i in range(0, 10):
    mcmc = pints.MCMCController(log_posterior, 3, xs,method=pints.HaarioBardenetACMC)
    mcmc.set_parallel(True)
    mcmc.set_max_iterations(30000)
    chains=mcmc.run()
    rhat_mean=np.mean(pints.rhat_all_params(chains))
    #pints.plot.trace(chains)
    #plt.show()
    if rhat_mean<1.2:
        Electrode="Red"
        Electrode_save="Red"
        folder="Noramp"
        Method="1"
        sim_options="1e-1M"
        run="MCMC_runs"
        save_file=filename=("_").join([folder,Method, sim_options, "MCMC"])+".run3"
        filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run])
        f=open(filepath+"/"+filename, "wb")
        np.save(f, chains)
        f.close()
        break
