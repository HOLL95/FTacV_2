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


file_list=["Noramp__1e0M_1_cv.run1", "Noramp__1e0M_2_cv.run1","Noramp__1e0M_3_cv.run1"]
param_list=[
            [0.4437541435808381, 0.04672244339040179, 169.82717561197694, 1499.9999999884658, 3.326734620433954e-05, -0.0024550352309892637, 0.003444105740388075, 1.6931489825942052e-10, 8.94186753976662, 4.218129657177546, 5.5734366297761, 0.5],
            [0.4480443313256469, 0.04578316505621488, 186.54661732132573, 1499.999998740905, 3.2749457208458154e-05, -0.0037426290649702418, 0.0033541225901685782, 1.5457625945960093e-10, 8.941954712163895, 4.2264148667766035, 5.54918407688431, 0.5],
            [0.44902174254951466, 0.04649482119604907, 167.77163066745334, 1499.9999899319125, 3.244942640679048e-05, -0.0015943375724840197, 0.003256602506039867, 1.5805338795893364e-10, 8.941715397492652, 4.227414643240183, 5.562298818419109, 0.5],
]
for i in range(2, len(file_list)):
    file=file_list[i]
    method="timeseries"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
    param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    param_vals=param_list[i]
    #noramp_results.dim_dict["v_nondim"]=True
    #param_vals[2]=150
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
    plt.plot(voltage_results, cmaes_time)
    plt.plot(voltage_results, current_results)
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
    noramp_results.param_bounds["Ru"]=[0, 2000]
    noramp_results.param_bounds["k_0"]=[0, 250]
    noramp_results.param_bounds["E0_mean"]=[0.3, 0.6]
    noramp_results.param_bounds["alpha"]=[0.3, 0.7]
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
    num_runs=10
    scores=np.ones(num_runs)*10
    for j in range(0, num_runs):
        current_min=min(scores)
        mcmc = pints.MCMCController(log_posterior, 3, xs,method=pints.HaarioBardenetACMC)
        mcmc.set_parallel(True)
        mcmc.set_max_iterations(30000)
        chains=mcmc.run()
        rhat_mean=np.mean(pints.rhat_all_params(chains[:, 20000:, :]))
        print(pints.rhat_all_params(chains[:, 20000:, :]))
        scores[j]=rhat_mean
        #k_rhat=pints.rhat_all_params(chains[:, 20000:, :])[2]
        #pints.plot.trace(chains)
        #plt.show()
        if rhat_mean<1.1:
            Electrode_save="Red"
            run2="MCMC_runs"
            save_file=file+"_MCMC_3a"
            filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run2])
            f=open(filepath+"/"+save_file, "wb")
            np.save(f, chains)
            f.close()
            break
