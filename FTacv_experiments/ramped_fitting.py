import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_ramped
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_ramped  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]

dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Black"
Method ="O_Method"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        voltages=np.loadtxt(path+"/"+data)
dec_amount=64
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=voltages[0::dec_amount, 1]

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
    'gamma': 1e-10,
    "original_gamma":1e-10,          # (surface coverage per unit area)
    'k_0': 100.0, #(reaction rate s-1)
    "E0_mean":0.2,
    "E0_std": 0.09,
    "k0_shape":0.954,
    "k0_loc":100,
    "k0_scale":50,
    "k0_range":1e3,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 50
}
likelihood_options=["timeseries", "fourier"]
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "test": False,
    "likelihood":likelihood_options[0],
    "dispersion":True,
    "dispersion_bins":16,
    "label": "cmaes",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,7,1),
    "experiment_time": time_results,
    "experiment_current": current_results,
    "experiment_voltage":voltage_results,
    "bounds_val":2000,
    "signal_length":len(time_results),
}


ramp_fit=single_electron(param_list, simulation_options, other_values)
time_results=ramp_fit.other_values["experiment_time"]
current_results=ramp_fit.other_values["experiment_current"]
voltage_results=ramp_fit.other_values["experiment_voltage"]
likelihood_func=ramp_fit.kaiser_filter(current_results)
ramp_fit.pass_extra_data(current_results, likelihood_func)
ramp_fit.optim_list=["gamma"]
ramp_fit.simulation_options["dispersion"]=False
normal=ramp_fit.test_vals([1e-10], likelihood="timeseries", test=False)
ramp_fit.simulation_options["dispersion"]=True
disped=ramp_fit.test_vals([1e-10], likelihood="timeseries", test=False)
plt.plot(normal)
plt.plot(disped)
plt.show()
ramp_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl','gamma', 'omega', 'phase', 'alpha']
param_bounds={
    'E_0':[0.1, 0.4],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 3e3],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-6], #(capacitance parameters)
    'CdlE1': [-1, 1],#0.000653657774506,
    'CdlE2': [0,0.1],#0.000245772700637,
    'CdlE3': [0,0.1],#1.10053945995e-06,
    'gamma': [1e-11,9e-10],
    'k_0': [0, 1e4], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    'k_0': 100.0, #(reaction rate s-1)
    "E0_mean":[0.0, 0.5],
    "E0_std": [0.01, 0.3],
    "k0_shape":[0,2],
    "k0_loc":[9e3, 1e4],
    "k0_scale":[0, 4e3],
    "k0_range":[1e2, 1e4],
    'phase' : [0, 2*math.pi]
}

ramp_fit.optim_list=['E0_mean', "E0_std",'k0_shape',"k0_loc", "k0_scale", "k0_range", 'Ru', 'Cdl',"CdlE1",'gamma', 'omega', 'phase', 'alpha']
param_boundaries=np.zeros((2, ramp_fit.n_parameters()))
for i in range(0, ramp_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[ramp_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[ramp_fit.optim_list[i]][1]

ramp_fit.define_boundaries(param_boundaries)
harm_class=harmonics(other_values["harmonic_range"], ramp_fit.nd_param.omega*ramp_fit.nd_param.c_T0, 0.1)
data_harmonics=harm_class.generate_harmonics(time_results, current_results)
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics,"abs", "numerical", "data")
#harm_class.harmonics_and_time(time_results, test_harmonics, test, data_time=current_results, harmonics2=data_harmonics,label1="numerical", label2="data", title="Black")
likelihood_func=ramp_fit.kaiser_filter(current_results)
cmaes_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
#cmaes_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, fourier_test)
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
x0=abs(np.random.rand(ramp_fit.n_parameters()))#ramp_fit.change_norm_group(means[1:], "norm")#abs(np.random.rand(ramp_fit.n_parameters()))

for i in range(0, 1):
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
#    x0=found_parameters
#    print found_parameters
print folder

cmaes_results=ramp_fit.change_norm_group(found_parameters, "un_norm")
cmaes_time=ramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
plt.plot(cmaes_time)
plt.plot(current_results)
plt.show()
print list(cmaes_results)
"""
cmaes_time_prediction=ramp_fit.simulate(found_parameters,frequencies, "optimise", "timeseries", "yes" )
test_harmonics=harm_class.generate_harmonics(time_results, cmaes_time_prediction)
harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics,"abs", "numerical", "data")
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)
error=np.std(np.subtract(cmaes_prediction, likelihood_func))
#error=np.std(np.subtract(cmaes_time_prediction, current_results))
mcmc_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, likelihood_func)
#mcmc_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
updated_lb=np.append(cmaes_results*0.5, [0])#found_parameters[3]*0.97,
updated_ub=np.append(cmaes_results*1.5, [max(likelihood_func)])#found_parameters[3]*1.03,
updated_boundaries=[updated_lb, updated_ub]
updated_boundaries=np.sort(updated_boundaries, 0)
ramp_fit.define_boundaries(updated_boundaries)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                updated_boundaries[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(cmaes_results, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
ramp_fit.label="MCMC"
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()
pints.plot.trace(chains)
plt.show()
filename="ramped_results_2"
f=open(filename, "w")
np.save(f, chains)
f.close()
"""
