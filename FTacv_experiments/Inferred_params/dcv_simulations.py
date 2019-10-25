import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy
import copy
import pints
import pints.plot
dir_path = os.path.dirname(os.path.realpath(__file__))
slash_idx=[i for i in range(len(dir_path)) if dir_path[i]=="/"]
one_above=dir_path[:slash_idx[-1]]
sys.path.insert(1, one_above)
from single_e_class_unified import single_electron
from harmonics_plotter import harmonics
Electrode="Yellow"
Run="Run_2"
path=("/").join([dir_path , Electrode, Run])
files=os.listdir(path)#
file_numbers=[]
counter=1
experiment=str(1)
repeat=str(1)
file1="Noramp_1_cv_low_ru.pkl"


def find(name, path, Electrode, skiprows=0):
    for root, dirs, files in os.walk(path):
        if Electrode in dirs:
            files=os.listdir(root+"/"+Electrode)
            if name in files:
                print(name)
                return np.loadtxt(root+"/"+Electrode+"/"+name, skiprows=skiprows)
param_dict={}
filename_1=["Ramped_3_cv_high_ru.ts"]
param_results=[]
labels=["Ramped"]
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"]
optim_dict={}



for names in filename_1:
    result=single_electron(path+"/"+names)
    param_results.append([result.save_dict["params"][0][result.save_dict["optim_list"].index(key)] if  (key in result.save_dict["optim_list"]) else result.dim_dict[key] for key in master_optim_list])
exp_currents=[]
blank_currents=[]
exp_names=[]
for i in range(1, 4):
    for j in range(1, 4):
        #plt.subplot(2,3,i)
        dcv_results=("_").join([Electrode, "Electrode", "Direct_CV", str(i), str(j)])
        dcv_blank=("_").join(["Blank", Electrode, "Electrode", "Direct_CV", str(i), str(j)])
        exp_results=find(dcv_results, one_above, Electrode,1)
        blank_results=find(dcv_blank, one_above, Electrode,1)
        blank_currents.append(exp_results[:,2])
        exp_currents.append(exp_results[:,2]-blank_results[:,2])
        exp_names.append(("_").join([Electrode, "DCV", str(i), str(j)]))
voltages=exp_results[:,1]
plt.plot(exp_currents[4])
plt.show()
dcv_params={
    "E_0":0.25,
    'E_start': -0.15, #(starting dc voltage - V)
    'E_reverse': 0.6,    #  (reverse dc voltage - V)
    'omega':8.94,  #    (frequency Hz)
    'd_E': 0,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 30e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 0, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,
    "original_gamma":1e-10,          # (surface coverage per unit area)
    'k_0': 1.0, #(reaction rate s-1)
    "E0_mean":0.2,
    "E0_std": 0.09,
    "k0_shape":0.954,
    "k0_loc":100,
    "k0_scale":50,
    "k0_range":1e3,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    "cap_phase":0,
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 600
}

simulation_options_dcv={
    "no_transient":5,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":20,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"dcv",
    "label": "MCMC",
    "optim_list":[]
}

dcv_other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(2,6,1)),
    "experiment_time": exp_results[:,0],
    "experiment_current":exp_currents[4], #noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":exp_results[:,1],#noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":len(exp_results[:,1]),
}
dcv_results=single_electron(None, dcv_params, simulation_options_dcv, dcv_other_values, result.param_bounds)
dcv_results.def_optim_list(["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma",'omega',"cap_phase","phase","alpha"])
time_results=dcv_results.other_values["experiment_time"]
current_results=dcv_results.other_values["experiment_current"]
voltage_results=dcv_results.other_values["experiment_voltage"]
voltages=voltage_results*dcv_results.nd_param.c_E0
results_dict={}

for i in range(3, 4):
    #current_results=exp_currents[i]

    cdl_idx=tuple(np.where((voltages<0.15) | (voltages>0.35)))
    farad_idx=tuple(np.where((voltages>0.15) & (voltages<0.35)))
    cdl_only=copy.deepcopy(current_results)
    times=time_results
    cdl_times=times#[cdl_idx]
    cdl_current=cdl_only#[cdl_idx]
    cdl_voltage=voltages[cdl_idx]
    halfway=np.where(voltages==max(voltages))[0][0]
    counter=1
    plt.subplot(2,4,4)
    plt.plot(times, current_results, alpha=0.5, label="data")
    plt.subplot(2,4,8)
    plt.plot(times, current_results, alpha=0.5, label="data")
    for q in np.arange(5, 3, -1):
        for k in range(1, 3):
            if k==1:
                cdl_times=times#[cdl_idx]
                cdl_current=cdl_only#[cdl_idx]
                n=0
            else:
                cdl_times=times[cdl_idx]
                cdl_current=cdl_only[cdl_idx]
                n=4
            interped=np.interp(times, cdl_times, cdl_current)
            interp1=interped[:halfway]
            interp2=interped[halfway:]
            time2=time_results[halfway:]
            time1=time_results[:halfway]
            fitted1=np.poly1d(np.polyfit(time1, interp1, q))
            fitted2=np.poly1d(np.polyfit(time2, interp2, q))
            plt.subplot(2,4,counter+n)
            total=np.append(fitted1(time1),fitted2(time2))
            plt.plot(voltages, np.subtract(current_results, total))
            plt.title("Order "+str(q))
            plt.subplot(2,4,4*k)
            plt.plot(times,  np.append(fitted1(time1),fitted2(time2)), label="Order "+str(q))
            plt.legend()
            plt.xlabel("Nondim time")
            plt.ylabel("Nondim current")
        counter+=1

plt.show()

print(dcv_results.param_bounds)
fitted_blank=np.append(fitted1(time1), fitted2(time2))
objective_func=np.subtract(current_results, fitted_blank)


#dcv_results=single_electron(None, dcv_params, simulation_options_dcv, dcv_other_values, result.param_bounds)
dcv_results.param_bounds["k_0"]=[0, 1e4]
dcv_results.param_bounds["Ru"]=[0, 1e4]
dcv_results.def_optim_list(["E0_mean", "E0_std","k_0","gamma","Ru","alpha"])#, "Cdl", "CdlE1", "CdlE2"])
dcv_results.simulation_options["likelihood"]="timeseries"
dcv_results.simulation_options["label"]="cmaes"
dcv_results.simulation_options["test"]=False

cmaes_problem=pints.SingleOutputProblem(dcv_results, dcv_results.other_values["experiment_time"], objective_func)
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(dcv_results.optim_list))], [np.ones(len(dcv_results.optim_list))])
print([np.zeros(len(dcv_results.optim_list))], [np.ones(len(dcv_results.optim_list))])
x0=abs(np.random.rand(dcv_results.n_parameters()))#dcv_results.change_norm_group(gc4_3_low_ru, "norm")
cmaes_fitting=pints.OptimisationController(score, x0, sigma0=None, boundaries=CMAES_boundaries, method=pints.CMAES)
cmaes_fitting.set_max_unchanged_iterations(iterations=200, threshold=1e-3)
cmaes_fitting.set_parallel(True)
found_parameters, found_value=cmaes_fitting.run()
cmaes_results=dcv_results.change_norm_group(found_parameters, "un_norm")
print(list(found_parameters))
print(list(cmaes_results))
cmaes_time=dcv_results.test_vals(cmaes_results, "timeseries")
plt.plot(dcv_results.other_values["experiment_time"], dcv_results.other_values["experiment_current"], label="data")
plt.plot(dcv_results.other_values["experiment_time"], objective_func, label="objective function")
plt.plot(dcv_results.other_values["experiment_time"], cmaes_time, label="simulations")
plt.show()
dcv_results.simulation_options["label"]="MCMC"
mcmc_problem=cmaes_problem=pints.SingleOutputProblem(dcv_results, dcv_results.other_values["experiment_time"], objective_func)

error=np.std(np.abs(np.subtract(cmaes_time, objective_func)))
updated_lb=np.append([0.2, 0, 0, 9e-12, 0, 0.3],0*error)
updated_ub=np.append([0.3, 0.1, 1e4, 2e-10, 1e4, 0.7], 10*error)
updated_b=[updated_lb, updated_ub]
updated_b=np.sort(updated_b, axis=0)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_b[0], updated_b[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(cmaes_results, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_parallel(True)
mcmc.set_max_iterations(15000)
chains=mcmc.run()
pints.plot.trace(chains)
plt.show()
inferred_params=np.zeros(len(dcv_results.optim_list))
for i in range(0, len(dcv_results.optim_list)):
    inferred_params[i]=np.mean(chains[:, 5000:, i])
mcmc_results=dcv_results.test_vals(inferred_params, "timeseries")
print(list(inferred_params))
plt.plot(times, mcmc_results)
plt.plot(times, objective_func)
plt.show()
save_file="DCV_MCMC_2_1"
f=open(("/").join([dir_path,Electrode,"DCV", "MCMC_runs", save_file]), "wb")
np.save(f, chains)
f.close()
