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
dir_path = os.path.dirname(os.path.realpath(__file__))
types=["current", "voltage"]
exp="Experimental-221119"
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
Electrode="Red"
folder="Noramp"
concs=["_1e0M_"]

for lcv_2 in range(0, len(concs)):
    total_voltages=[]
    total_time=[0]
    total_current=[]
    for lcv_1 in range(1, 4):
        Method =concs[lcv_2]+str(lcv_1)+"_cv"
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
        kernel_size=8
        #current_results=np.convolve(results[:, 1], np.ones((kernel_size,))/kernel_size, mode="valid")
        #current_results1=current_results[0::dec_amount]

        current_results1=results[0::dec_amount, 1]
        print(len(current_results1))
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
            'sampling_freq' : (1.0/200),
            'phase' : 3*(math.pi/2),
            "time_end": None,
            'num_peaks': 20
        }
        solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
        likelihood_options=["timeseries", "fourier"]
        time_start=0.1/(param_list["omega"])
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
            'k_0': [50, 1e3], #(reaction rate s-1)
            'alpha': [0.3, 0.7],
            "cap_phase":[0, 2*math.pi],
            "E0_mean":[0.15, 0.25],
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
        total_voltages=np.append(total_voltages, voltage_results)
        total_current=np.append(total_current, current_results)
        total_time=np.append(total_time, np.add(time_results, total_time[-1]))
    total_time=np.delete(total_time, 0)
    other_values["experiment_time"]=total_time
    other_values["experiment_voltage"]=total_voltages
    other_values["experiment_current"]=total_current

    noramp_fit.nd_param.time_end=total_time[-1]
    noramp_fit.times()
    noramp_fit.simulation_options["interpolant"]=total_time
    noramp_fit.simulation_options["no_transient"]=False
    #noramp_fit.time_idx=0
    #noramp_fit.time_idx=0
    noramp_fit.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase","alpha"])
    true_data=total_current
    inf_params=[0.15000000018944165, 0.10585753483158893, 58.647835768143864, 3086.2298518853395, 5.642632260908315e-05, 0.05396499954590106, -0.00039767166962080894, 3.0365465573369006e-10, 8.944854891819027, 4.248411764775375, 0.18772534522877896, 0.5]

    test=noramp_fit.test_vals(inf_params, "timeseries")
    test=np.interp(total_time, noramp_fit.time_vec[noramp_fit.time_idx:], test)

    plt.plot(total_time, total_current, label="0.1M NaCl 1-3 concatenated")
    plt.plot(total_time, test, label="Sim")
    plt.xlabel("Nondim time")
    plt.ylabel("Nondim current")
    plt.legend()
    plt.show()
    #plt.plot(voltage_results, current_results)
    #plt.show()
    fourier_arg=noramp_fit.kaiser_filter(current_results)
    if simulation_options["likelihood"]=="timeseries":
        cmaes_problem=pints.SingleOutputProblem(noramp_fit, total_time, true_data)
    elif simulation_options["likelihood"]=="fourier":
        dummy_times=np.linspace(0, 1, len(fourier_arg))
        cmaes_problem=pints.SingleOutputProblem(noramp_fit, dummy_times, fourier_arg)
    score = pints.SumOfSquaresError(cmaes_problem)#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
    CMAES_boundaries=pints.RectangularBoundaries(list([np.zeros(len(noramp_fit.optim_list))]), list([np.ones(len(noramp_fit.optim_list))]))
    noramp_fit.simulation_options["label"]="cmaes"
    #noramp_fit.simulation_options["test"]=True
    num_runs=20
    param_mat=np.zeros((num_runs,len(noramp_fit.optim_list)))
    score_vec=np.zeros(num_runs)
    for i in range(0, num_runs):
        x0=abs(np.random.rand(noramp_fit.n_parameters()))#noramp_fit.change_norm_group(gc4_3_low_ru, "norm")
        print(noramp_fit.change_norm_group(x0, "un_norm"))
        cmaes_fitting=pints.OptimisationController(score, x0, sigma0=None, boundaries=CMAES_boundaries, method=pints.CMAES)
        cmaes_fitting.set_max_unchanged_iterations(iterations=200, threshold=1e-3)
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
        #plt.plot(total_time, cmaes_time)
        #plt.plot(total_time, total_current)
        #plt.show()
        cmaes_fourier=noramp_fit.test_vals(cmaes_results, likelihood="fourier", test=False)
        param_mat[i,:]=cmaes_results
        score_vec[i]=found_value

    idx=[i[0] for i in sorted(enumerate(score_vec), key=lambda y:y[1])]
    save_params=param_mat[idx[0:3], :]
    best_idx=np.where(score_vec==min(score_vec))
    best_idx=best_idx[0][0]
    cmaes_results=param_mat[best_idx,:]
    cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
    Electrode_save=extra+Electrode
    run="Run_1"
    if "k0_shape" in noramp_fit.optim_list:
        sim_options=resistance_type+"_"+"k0_disp"
    else:
        sim_options=resistance_type
    filename=("_").join([folder,Method])+".concatenated_data"
    filepath=("/").join([dir_path, "Inferred_params", Electrode_save, run])
    noramp_fit.save_state(results_dict, filepath, filename, save_params)
