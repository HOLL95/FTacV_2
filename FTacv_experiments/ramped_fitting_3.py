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
exp="Experimental-120919"
bla="Blank-110919"
resistances=["high_ru", "low_ru", "fixed_ru"]
ru_upper_bound=[1e4, 85, 50]
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
folder="Ramped"

Method ="2_cv"
type="current"
type2="voltage"
path=("/").join([dir_path, data_path, folder, Electrode])
files= os.listdir(path)
print(files)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
        print(data)
    elif (Method in data)  and (type2 in data):
        voltages=np.loadtxt(path+"/"+data)

dec_amount=1
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
current_results1=results[0::dec_amount, 1]
time_results1=results[0::dec_amount, 0]
voltage_results1=voltages[0::dec_amount, 1]
results_dict={"experiment_voltage": voltage_results1,
                "experiment_time": time_results1,
                "experiment_current": current_results1,}
plt.plot(time_results1, current_results1)
plt.show()
param_list={
    "E_0":0.2,
    'E_start': -180e-3, #(starting dc voltage - V)
    'E_reverse': 620e-3,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 29.8e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 20.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-6, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    "original_gamma": 1e-10,
    'gamma': 1e-10,
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
    'phase' : 0,
    "time_end": None,
    'num_peaks': 50,
    "v_nondim":True
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
time_start=5/(param_list["omega"])
start_idx=np.where(time_results1>time_start)
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":5,
    "test": False,
    "GH_quadrature":True,
    "phase_only":True,
    "method": "ramped",
    "likelihood":likelihood_options[1],
    "numerical_method": solver_list[1],
    "label": "MCMC",
    "optim_list":[]
}
for j in range(3, 7, 1):
    other_values={
        "filter_val": 0.1,
        "harmonic_range":list(range(4,7,1)),
        "experiment_time": time_results1,
        "experiment_current": current_results1,
        "experiment_voltage":voltage_results1,
        "bounds_val":20,
        "signal_length":int(len(time_results1))
    }
    param_bounds={
        'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
        'omega':[0.95*param_list['omega'],1.05*param_list['omega']],#8.88480830076,  #    (frequency Hz)
        'Ru': [0, 3000],  #     (uncompensated resistance ohms)
        'Cdl': [0,1e-4], #(capacitance parameters)
        'CdlE1': [-0.05,0.15],#0.000653657774506,
        'CdlE2': [-0.01,0.01],#0.000245772700637,
        'CdlE3': [-0.01,0.01],#1.10053945995e-06,
        'gamma': [1e-11,1e-9],
        'k_0': [1, 200], #(reaction rate s-1)
        'alpha': [0.4, 0.6],
        "cap_phase":[0, 2*math.pi],
        "E0_mean":[0.2, 0.3],
        "E0_std": [0.001, 0.2],
        "k0_shape":[0,2],
        "k0_loc":[0, 1e3],
        "k0_scale":[0,1e3],
        "k0_range":[1e2, 1e4],
        'phase' : [0, 2*math.pi]
    }
    ramp_fit=single_electron(None, param_list, simulation_options, other_values, param_bounds)
    sim_volts=ramp_fit.define_voltages()
    time_results=ramp_fit.other_values["experiment_time"]
    current_results=ramp_fit.other_values["experiment_current"]
    voltage_results=ramp_fit.other_values["experiment_voltage"]
    print(len(current_results))
    ramp_harm=harmonics(other_values["harmonic_range"], ramp_fit.nd_param.omega*ramp_fit.nd_param.c_T0, 0.08)
    data_harms=ramp_harm.generate_harmonics(time_results, current_results)
    #ramp_fit.def_optim_list(["E0_mean","E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"])
    ramp_fit.def_optim_list(["E0_mean", "E0_std","k_0" ,"Ru","Cdl","CdlE1", "CdlE2","gamma",'omega'])
    #ramp_fit.dim_dict["gamma"]=0
    noramp_params=[0.23363378517047495, 0.03481010462713212, 125.2418680556816, 630.0772611187369, 7.694171200331618e-05, 0.003209611999861764, -0.0004263185805571494, 7.476933579891946e-11, 8.940473342512918]
    ramp_fit.def_optim_list(["E0_mean", "E0_std","k_0" ,"Ru","Cdl","CdlE1", "CdlE2","gamma",'omega', "alpha"])
    inferred_params=[0.21363378517047495, 0.0426544633855209, 122.16602260400302, 742.0285824092344, 8.212777750912338e-05, 0.00295011584929738, -0.0005139373308230138, 7.130227126430568e-11, 8.875659818458482, 0.4274370247643113, 0.5146469208293426]






    print(ramp_fit.param_bounds)
    #ramp_fit.def_optim_list(["Ru","Cdl","CdlE1", "CdlE2",'omega',"phase","cap_phase"])
    #ramp_fit.dim_dict["gamma"]=0
    """
    if ru_pick==2:
        ramp_fit.def_optim_list(["E0_mean", "E0_std","k_0" ,"Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase", "alpha"])
        ramp_fit.dim_dict["Ru"]=65
    """
        #ramp_fit.dim_dict["alpha"]=0.5
    true_data=current_results
    #ramp_fit.simulation_options["fourier_scaling"]=False
    if simulation_options["likelihood"]=="timeseries":
        cmaes_problem=pints.SingleOutputProblem(ramp_fit, time_results, true_data)
    elif simulation_options["likelihood"]=="fourier":
        fourier_arg=ramp_fit.kaiser_filter(current_results)
        ramp_fit.secret_data_fourier=fourier_arg
        dummy_times=np.linspace(0, 1, len(fourier_arg))
        cmaes_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, fourier_arg)
    score = pints.SumOfSquaresError(cmaes_problem)#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
    CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(ramp_fit.optim_list))], [np.ones(len(ramp_fit.optim_list))])
    ramp_fit.simulation_options["label"]="cmaes"
    ramp_fit.simulation_options["test"]=False
    num_runs=5
    param_mat=np.zeros((num_runs,len(ramp_fit.optim_list)))
    score_vec=np.zeros(num_runs)
    for i in range(0, num_runs):
        noramp_params=[0.23863378517047495, 0.04481010462713212, 125.2418680556816, 512.0772611187369, 7.694171200331618e-05, 0,0, 7.176933579891946e-11, 8.840473342512918,0, 0.5999996422725197]
        noramp_params=[0.22774874566415484, 0.046075057493608136, 123.33007397100599, 873.5412252656006, 9.999999955830733e-05, -0.04999999026295966, 0.0028305800048256807, 7.378704532562238e-11, 8.88472746776295, 0.5900000005069297]
        noramp_params=[0.22925918516708654, 0.04595696579195954, 123.33007397100599, 873.5412252656006, 3.3412012933121965e-05, 0.057928207116134806, -0.0021217096115628917, 7.178042062464878e-11, 8.884751771027027, 0.43751189633466997]
        x0=ramp_fit.change_norm_group(noramp_params, "norm")#abs(np.random.rand(ramp_fit.n_parameters()))#ramp_fit.change_norm_group(gc4_3_low_ru, "norm")
        cmaes_fitting=pints.OptimisationController(score, x0, sigma0=None, boundaries=CMAES_boundaries, method=pints.CMAES)
        cmaes_fitting.set_max_unchanged_iterations(iterations=200, threshold=1e-3)
        if "E0_mean" in ramp_fit.optim_list and "k0_loc" in ramp_fit.optim_list:
            cmaes_fitting.set_parallel(False)
        else:
            cmaes_fitting.set_parallel(True)
        found_parameters, found_value=cmaes_fitting.run()
        found_value=1000
        cmaes_results=ramp_fit.change_norm_group(found_parameters, "un_norm")
        #cmaes_results=noramp_params
        print(list(cmaes_results))
        cmaes_time=ramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
        cmaes_fourier=ramp_fit.test_vals(cmaes_results, likelihood="fourier", test=False)
        param_mat[i,:]=cmaes_results
        score_vec[i]=found_value
        #plt.plot(time_results, current_results, alpha=0.6)
        #plt.plot(time_results, cmaes_time, alpha=0.6)
        #plt.show()
        sim_harms=ramp_harm.generate_harmonics(time_results, cmaes_time)
        #for j in range(0, len(data_harms)):
    #        plt.subplot(len(data_harms), 1, j+1)
    #        plt.plot(time_results, (data_harms[j,:]))
    #        plt.plot(ramp_fit.time_vec, (sim_harms[j,:]), alpha=0.5)
    #    plt.show()
        #plt.plot(time_results, cmaes_time)
        #plt.plot(time_results, true_data, alpha=0.5)
        #plt.plot(time_results, np.subtract(true_data, cmaes_time))
        #plt.show()
idx=[i[0] for i in sorted(enumerate(score_vec), key=lambda y:y[1])]
save_params=param_mat[idx[0:3], :]
Electrode_save=extra+Electrode
if "k0_shape" in ramp_fit.optim_list:
    sim_options=resistance_type+"_"+"k0_disp"
else:
    sim_options=resistance_type
filename=("_").join([folder,Method, sim_options])+".blank"
filepath=("/").join([dir_path, "Inferred_params", Electrode_save])
ramp_fit.save_state(results_dict, filepath, filename, save_params)
best_idx=np.where(score_vec==min(score_vec))
best_idx=best_idx[0][0]
cmaes_results=param_mat[best_idx,:]
print(list(cmaes_results))
cmaes_time=ramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
