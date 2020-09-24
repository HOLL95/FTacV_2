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
run="Run_7"
def RMSE(series1, series2):
    return np.sqrt((np.sum(1/(len(series1))*np.power(np.subtract(series1, series2),2))))
final_osc=[25, 20, 15, 10]
other_files=["6", "9"]
param_matrix=[[0.24696678941258246, 0.02150713711591799, 99.45832461312207, 774.704213209354, 7.701956070978554e-05, 0.003017010623506844, -0.00039055262925760496, 7.758858523468079e-11, 8.940640920538685, 4.410033377066397,5.1069657955459,  0.5825105457091572, 0.14328732063771737] ,
        [0.2438679560342186, 0.04933431717529572, 191.31132541897506, 497.5055066560729, 7.69338292657605e-05, 0.0031326600706418906, -0.0004312564011993529, 7.271953365051911e-11, 8.94052385943315, 4.3351579045955795,4.925339460865147,  0.6496639586748201, 0.17177231234727677] ,
        [0.24430949230426569, 0.04731795471469772, 176.34896534285173, 521.4259532758771, 7.647752507363084e-05, 0.002994039070869009, -0.0004190299720139099, 7.275919967176139e-11, 8.940534898414818, 4.342816983799178, 4.943663135736721, 0.6362752123549957, 0.1720460446880876] ,
        [0.24390113149368509, 0.0457048751819453, 161.20949406256256, 536.934970450984, 7.621126298448553e-05, 0.0028379095150146416, -0.000413526916786134, 7.281036191448345e-11, 8.940539372196708, 4.34670259354896, 4.957545393990694, 0.6227137060882207, 0.1713086829766499] ,
        [0.24264240746946078, 0.04862334843813994, 174.5319444179432, 501.7381191076415, 7.621609295572e-05, 0.0025967801809421154, -0.00041475793632247265, 7.236585251042315e-11, 8.94052009496047, 4.336548651366052, 4.934251569752124, 0.6291348730960813, 0.17453584870953875] ,
        [0.24160781666577544, 0.05036058151555205, 183.91624120363147, 478.0527573517239, 7.62482655530355e-05, 0.002429757046640779, -0.00041668758601426555, 7.209545158956022e-11, 8.940507420558243, 4.330200828730097, 4.9203855923331705, 0.633836591099598, 0.17606178914368914] ,
        [0.24122209427368013, 0.049981782635136945, 179.81828332106306, 480.1074385658014, 7.615225208615455e-05, 0.002305928382454525, -0.00041429576393527123, 7.19385317237139e-11, 8.940492881836837, 4.331293568942899, 4.923830104418481, 0.6292049920997411, 0.17637075996545973] ,
        [0.24100346908927325, 0.05022099349663932, 178.28549652948317, 477.9310528442049, 7.59455358542439e-05, 0.002273221478669929, -0.00041353404857902287, 7.160330874703819e-11, 8.940504038013483, 4.329624942864246, 4.92163560153165, 0.6277651371395521, 0.17665475180369178] ,
        [0.2412440781552388, 0.050226126358631266, 180.12768848217763, 485.192042495465, 7.548546629894573e-05, 0.0021505398969171585, -0.00040389051055751864, 7.11413102722647e-11, 8.940519733026564,  4.3302796214687795,4.9224469580118715, 0.6290049217531675, 0.17737634083075238],
        [0.24055747088306997, 0.05047523191357855, 176.56116945811587, 474.43867201525205, 7.584083348640005e-05, 0.0021028609527935505, -0.0004094927674505519, 7.146018738751102e-11, 8.940505868813535,  4.327745982059323, 4.920060562270928,0.6246797783009037, 0.1779493870477656]]
for i in (range(1, 11)):
    #file="Noramp_"+str(i)+"_cv_high_ru_alpha_disp"

    file="Noramp_"+str(i)+"_cv_high_ru_alpha_disp"

    #file="Noramp_2_cv_high_ru_alpha_0.61"
    method="timeseries"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase", "phase","alpha_std", "alpha_mean"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
    param_vals=([noramp_results.save_dict["params"][1][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    noramp_results.def_optim_list(master_optim_list)
    ts=noramp_results.test_vals(param_vals, "timeseries")
    err=RMSE(ts, noramp_results.other_values["experiment_current"])
    print_vals=param_vals+[noramp_results.i_nondim(err)*1e6]
    print(print_vals, ",")
   # param_vals=param_matrix[i-1]
   # print(noramp_results.dim_dict["sampling_freq"], noramp_results.time_vec[1]-noramp_results.time_vec[0])
    #param_vals[0]=noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index("E0_mean")]
    #noramp_results.dim_dict["v_nondim"]=True
    #param_vals[2]=150
    #noramp_results.simulation_options["label"]="MCMC"
   # print(param_vals)
    noramp_results.dim_dict["alpha_mean"]=noramp_results.save_dict["params"][2][noramp_results.save_dict["optim_list"].index("alpha_mean")]
    noramp_results.simulation_options["GH_quadrature"]=True
   # print(noramp_results.simulation_options["dispersion_bins"])
   # print(noramp_results.save_dict["optim_list"])
   # noramp_results.simulation_options["dispersion_bins"]=[8,8]
   # del param_vals[-2]
    noramp_results.def_optim_list(master_optim_list)


    #noramp_results.simulation_options["likelihood"]="fourier"


    #plt.plot(noramp_results.other_values["experiment_time"], noramp_results.other_values["experiment_current"])
    #current_results=np.convolve(noramp_results.other_values["experiment_current"], np.ones((kernel_size,))/kernel_size, mode="valid")
    current_results=noramp_results.other_values["experiment_current"]
    #current_results=np.append(current_results,noramp_results.other_values["experiment_current"][-1] )
    voltage_results=noramp_results.other_values["experiment_voltage"]
    time_results=noramp_results.other_values["experiment_time"]

    #cmaes_time=noramp_results.test_vals(param_vals, method)

    method="timeseries"
    if method=="timeseries":
        fit_data=current_results
        fit_times=time_results
    elif method=="fourier":
        fit_data=noramp_results.kaiser_filter(current_results)
        fit_times=np.linspace(0, 1, len(fit_data))
    cmaes_time=noramp_results.test_vals(param_vals, "timeseries")
    error=RMSE(cmaes_time, current_results)
    param_vs=param_vals+[error]
    print(param_vs, ",")
    """
    print("Error", error, noramp_results.i_nondim(error)*1e6)
    mcmc_problem=pints.SingleOutputProblem(noramp_results, fit_times, fit_data)

    #updated_lb=np.append([noramp_results.param_bounds[key][0] for key in master_optim_list],0.75*error)
    #updated_ub=np.append([noramp_results.param_bounds[key][1] for key in master_optim_list], 1.25*error)


    noramp_results.param_bounds["Ru"]=[0, 1000]
    noramp_results.param_bounds["k_0"]=[0, 500]
    noramp_results.param_bounds["alpha_std"]=[1e-5, 0.3]
    #noramp_results.simulation_options["alpha_dispersion"]="normal"
    #noramp_results.dim_dict["alpha_mean"]=0.6
    noramp_results.dim_dict["alpha_std"]=1e-2
    noramp_results.def_optim_list(master_optim_list)
    updated_lb=np.append([noramp_results.param_bounds[key][0] for key in master_optim_list], 0.01*error)
    updated_ub=np.append([noramp_results.param_bounds[key][1] for key in master_optim_list], 10*error)
    mcmc_problem=pints.SingleOutputProblem(noramp_results, time_results, current_results)
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
        #alpha_index=noramp_results.optim_list.index("alpha")
        alpha_chain=[]
        mcmc.set_parallel(True)
        mcmc.set_max_iterations(100000)
        chains=mcmc.run()
        rhat_mean=np.mean(pints.rhat_all_params(chains[:, 50000:, :]))
        #for q in range(0, 2):
        #    alpha_chain=np.append(alpha_chain, chains[q, 30000:, alpha_index])
        Electrode_save="Yellow"
        run2="MCMC_runs/omega_nondim"
        save_file=file+"_MCMC_run_25"
        filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run2])
        #print(pints.rhat_all_params(chains[:, 20000:, :]))
        #k_rhat=pints.rhat_all_params(chains[:, 20000:, :])[2]
        #pints.plot.trace(chains)
        #plt.show()
        if rhat_mean<1.08:
            print("Good save")
            f=open(filepath+"/"+save_file, "wb")
            np.save(f, chains)
            f.close()
            break
        elif rhat_mean<min(scores):
            print("bad save")
            f=open(filepath+"/"+save_file, "wb")
            np.save(f, chains)
            f.close()
        scores[j]=rhat_mean

"""
