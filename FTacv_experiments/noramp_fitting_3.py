import numpy as np
import matplotlib.pyplot as plt
import math
import time
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
import os
import pickle
import pints
import sys
import select
dir_path = os.path.dirname(os.path.realpath(__file__))
types=["current", "voltage"]
exp="Experimental-120919"
bla="Blank-110919"
resistances=["high_ru", "low_ru", "fixed_ru"]
ru_upper_bound=[700, 1e4, 50]
ru_pick=0
resistance_type=resistances[ru_pick]
print(resistance_type)
exp_type=exp
if exp_type==bla:
    extra="Blank/"
else:
    extra=""
data_path="/experiment_data_2/"+exp_type
Electrode="Yellow"
folder="Noramp"
for lcv_1 in range(1, 11):
    Method =str(lcv_1)+"_cv"
    type="current"
    type2="voltage"
    path=("/").join([dir_path, data_path, folder, Electrode])
    files= os.listdir(path)
    for data in files:
        if (Method in data)  and (type in data):
            results=np.loadtxt(path+"/"+data)
            print(data)
        elif (Method in data)  and (type2 in data):
            voltages=np.loadtxt(path+"/"+data)

    dec_amount=32
    de=300e-3
    estart=260e-3-de
    ereverse=estart+2*de
    current_results1=results[0::dec_amount, 1]
    time_results1=results[0::dec_amount, 0]
    voltage_results1=voltages[0::dec_amount, 1]
    results_dict={"experiment_voltage": voltage_results1,
                    "experiment_time": time_results1,
                    "experiment_current": current_results1,}
    param_list={
        "E_0":0.2,
        'E_start': estart, #(starting dc voltage - V)
        'E_reverse': ereverse,
        'omega':8.94,#8.88480830076,  #    (frequency Hz)
        "original_omega":8.94,
        'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
        'v': 10.36e-3,   #       (scan rate s^-1)
        'area': 0.07, #(electrode surface area cm^2)
        'Ru': 1.0,  #     (uncompensated resistance ohms)
        'Cdl': 1e-5, #(capacitance parameters)
        'CdlE1': 0,#0.000653657774506,
        'CdlE2': 0,#0.000245772700637,
        'CdlE3': 0,#1.10053945995e-06,
        'gamma': 1e-10,
        "original_gamma":1e-10,        # (surface coverage per unit area)
        'k_0': 10, #(reaction rate s-1)
        'alpha': 0.5,
        "E0_mean":0.2,
        "E0_std": 0.09,
        "k0_shape":0.954,
        "k0_loc":100,
        "k0_scale":50,
        "k0_range":1e3,
        "cap_phase":0,
        'sampling_freq' : (1.0/400),
        'phase' : 3*(math.pi/2),
        "time_end": None,
        'num_peaks': 30
    }
    solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
    likelihood_options=["timeseries", "fourier"]
    time_start=2/(param_list["omega"])
    simulation_options={
        "no_transient":time_start,
        "numerical_debugging": False,
        "experimental_fitting":True,
        "dispersion":False,
        "dispersion_bins":16,
        "test": False,
        "method": "sinusoidal",
        "phase_only":False,
        "likelihood":likelihood_options[0],
        "numerical_method": solver_list[1],
        "label": "MCMC",
        "optim_list":[]
    }
    other_values={
        "filter_val": 0.5,
        "harmonic_range":list(range(3,9,1)),
        "experiment_time": time_results1,
        "experiment_current": current_results1,
        "experiment_voltage":voltage_results1,
        "bounds_val":20,
        "signal_length":int(1e4)
    }
    param_bounds={
        'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
        'omega':[0.95*param_list['omega'],1.05*param_list['omega']],#8.88480830076,  #    (frequency Hz)
        'Ru': [0, 1e4],  #     (uncompensated resistance ohms)
        'Cdl': [0,1e-4], #(capacitance parameters)
        'CdlE1': [-0.05,0.15],#0.000653657774506,
        'CdlE2': [-0.01,0.01],#0.000245772700637,
        'CdlE3': [-0.01,0.01],#1.10053945995e-06,
        'gamma': [1e-11,1e-9],
        'k_0': [10, 1e3], #(reaction rate s-1)
        'alpha': [0.4, 0.6],
        "cap_phase":[0, 2*math.pi],
        "E0_mean":[0.19, 0.25],
        "E0_std": [0.001, 0.2],
        "k0_shape":[0,2],
        "k0_loc":[0, 1e3],
        "k0_scale":[0,1e3],
        "k0_range":[1e2, 1e4],
        'phase' : [0, 2*math.pi]
    }
    noramp_fit=single_electron(None, param_list, simulation_options, other_values, param_bounds)
    noramp_fit.define_boundaries(param_bounds)
    time_results=noramp_fit.other_values["experiment_time"]
    current_results=noramp_fit.other_values["experiment_current"]
    voltage_results=noramp_fit.other_values["experiment_voltage"]
    noramp_fit.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase","alpha"])
    #noramp_fit.def_optim_list(["Ru","Cdl","CdlE1", "CdlE2",'omega',"phase","cap_phase"])
    #noramp_fit.dim_dict["gamma"]=0
    true_data=current_results
    #plt.plot(voltage_results, current_results)
    #plt.show()
    fourier_arg=noramp_fit.kaiser_filter(current_results)
    if simulation_options["likelihood"]=="timeseries":
        cmaes_problem=pints.SingleOutputProblem(noramp_fit, time_results, true_data)
    elif simulation_options["likelihood"]=="fourier":
        dummy_times=np.linspace(0, 1, len(fourier_arg))
        cmaes_problem=pints.SingleOutputProblem(noramp_fit, dummy_times, fourier_arg)
    score = pints.SumOfSquaresError(cmaes_problem)#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
    CMAES_boundaries=pints.RectangularBoundaries(list([np.zeros(len(noramp_fit.optim_list))]), list([np.ones(len(noramp_fit.optim_list))]))
    noramp_fit.simulation_options["label"]="cmaes"
    #noramp_fit.simulation_options["test"]=True
    num_runs=10
    param_mat=np.zeros((num_runs,len(noramp_fit.optim_list)))
    score_vec=np.ones(num_runs)*1e6
    for i in range(0, num_runs):
        x0=abs(np.random.rand(noramp_fit.n_parameters()))#noramp_fit.change_norm_group(gc4_3_low_ru, "norm")
        print(noramp_fit.change_norm_group(x0, "un_norm"))
        cmaes_fitting=pints.OptimisationController(score, x0, sigma0=None, boundaries=CMAES_boundaries, method=pints.CMAES)
        cmaes_fitting.set_max_unchanged_iterations(iterations=200, threshold=1e-7)
        if "E0_mean" in noramp_fit.optim_list and "k0_loc" in noramp_fit.optim_list:
            cmaes_fitting.set_parallel(False)
        else:
            cmaes_fitting.set_parallel(True)
        found_parameters, found_value=cmaes_fitting.run()
        cmaes_results=noramp_fit.change_norm_group(found_parameters, "un_norm")
        print(list(cmaes_results))
        cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
        print("omega", noramp_fit.nd_param.nd_omega, noramp_fit.nd_param.omega)
        #plt.subplot(1,2,1)
        #plt.plot(voltage_results, cmaes_time)
        #plt.plot(voltage_results, current_results)
        #plt.subplot(1,2,2)
        #plt.plot(time_results, noramp_fit.define_voltages()[noramp_fit.time_idx:])
        #plt.plot(time_results, voltage_results)
        #plt.show()
        cmaes_fourier=noramp_fit.test_vals(cmaes_results, likelihood="fourier", test=False)
        param_mat[i,:]=cmaes_results
        score_vec[i]=found_value
        print("Finish?")
        i, o, e = select.select( [sys.stdin], [], [], 5)
        if len(i) != 0:
            break


    idx=[i[0] for i in sorted(enumerate(score_vec), key=lambda y:y[1])]
    save_params=param_mat[idx[0:3], :]
    best_idx=np.where(score_vec==min(score_vec))
    best_idx=best_idx[0][0]
    cmaes_results=param_mat[best_idx,:]
    cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
    Electrode_save=extra+Electrode
    run="Run_4"
    if "k0_shape" in noramp_fit.optim_list:
        sim_options=resistance_type+"_"+"k0_disp"
    else:
        sim_options=resistance_type
    filename=("_").join([folder,Method, sim_options])+".run5"
    filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run])
    noramp_fit.save_state(results_dict, filepath, filename, save_params)
    error=np.std(abs(np.subtract(cmaes_time, current_results)))
    mcmc_problem=pints.SingleOutputProblem(noramp_fit, time_results, current_results)
    updated_lb=np.append([noramp_fit.param_bounds[x][0] for x in noramp_fit.optim_list],0.01*error)
    updated_ub=np.append([noramp_fit.param_bounds[x][1] for x in noramp_fit.optim_list],10*error)
    updated_b=[updated_lb, updated_ub]
    updated_b=np.sort(updated_b, axis=0)
    log_liklihood=pints.GaussianLogLikelihood(mcmc_problem)
    #print(noramp_fit.n_parameters(), len(updated_b[0]))
    log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
    #print(log_prior.n_parameters(), log_liklihood.n_parameters())
    log_posterior=pints.LogPosterior(log_liklihood, log_prior)
    #[(noramp_fit.param_bounds[x][1]+noramp_fit.param_bounds[x][0])/2 for x in noramp_fit.optim_list ]
    mcmc_parameters=cmaes_results
    mcmc_parameters=np.append(mcmc_parameters, error)
    xs=[mcmc_parameters,
        mcmc_parameters,
        mcmc_parameters
        ]
    noramp_fit.simulation_options["label"]="MCMC"
    for i in range(0, 10):
        mcmc = pints.MCMCController(log_posterior, 3, xs,method=pints.HaarioBardenetACMC)
        mcmc.set_parallel(True)
        mcmc.set_max_iterations(40000)
        chains=mcmc.run()
        rhat_mean=np.mean(pints.rhat_all_params(chains[:, 30000:, :]))
        #pints.plot.trace(chains)
        #plt.show()
        if rhat_mean<1.1:
            run="MCMC_runs/omega_nondim"
            save_file=filename=("_").join([folder,Method, sim_options, "MCMC"])+".run6"
            filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run])
            f=open(filepath+"/"+filename, "wb")
            np.save(f, chains)
            f.close()
            break
