import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_ramped
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from single_e_class_ramped  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]

dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Carbon"
Method ="GC4_1_ramp"
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
    "cap_phase":0,
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
    "dispersion":False,
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
#plt.plot(normal)
#plt.plot(disped)
#plt.show()
ramp_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl','gamma', 'omega', 'phase', 'alpha']
param_bounds={
    'E_0':[0.1, 0.4],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 3e3],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-3], #(capacitance parameters)
    'CdlE1': [-5, 5],#0.000653657774506,
    'CdlE2': [-0.1,0.1],#0.000245772700637,
    'CdlE3': [0,0.1],#1.10053945995e-06,
    'gamma': [1e-11,1e-9],
    'k_0': [0, 1e4], #(reaction rate s-1)
    'alpha': [0.2, 0.8],
    'k_0': [0, 1e4], #(reaction rate s-1)
    "E0_mean":[0.0, 0.5],
    "E0_std": [0.01, 0.3],
    "k0_shape":[0,5],
    "k0_loc":[1, 1e4],
    "k0_scale":[0,2e3],
    "k0_range":[1e2, 1e4],
    'phase' : [0, 2*math.pi]
}
ramp_fit.optim_list=[]

harm_class=harmonics(other_values["harmonic_range"], ramp_fit.nd_param.omega*ramp_fit.nd_param.c_T0, 0.1)
ramp_fit.optim_list=["E0_mean", "E0_std","k_0", 'Ru', 'Cdl',"CdlE1","CdlE2",'gamma', 'omega', 'phase', 'alpha']


ramp_means_black=[0.20504878429957712, 0.07357984171424121, 126.17588866800264, 24.24736830211999, 9.999999999876686e-07, 0.11886527586092166, 0,1.3175840584818795e-10, 8.959496500290681, 6.2831853071795845, 0.599999999999731]
ramp_means_black_2=[0.21200070197828386, 0.06888959981526988, 133.96563492653507, 40.08177557223102, 3.226207450320691e-06, -0.021487125154184827, 0.0017931237883237632, 1.2669148148700962e-10, 8.959483328711777, 6.283185307173828, 0.7999999999999661]
ramp_means_carbon=[0.23192913053278295, 0.07597303082211063,133.999986524228, 20, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 2.8654939556021e-10, 8.959351751379364, 6.077678909557666, 0.7999999909807196]
ramp_free_means_black=[0.20504878429957712, 0.04692985835905884, 773.0039074468887, 1.1172494386860095e-06, 6.253912022948444e-06, 0.46590284463560927, -0.020672008906663236, 9.665438282298918e-11, 8.94055300929529,6.2831853071795845,  0.10000008017506497]
ramp_free_means_black_2=[0.22026089333976873, 0.04776183826475387, 1226.0003897156193, 2.6091841820962103e-10, 6.311657164346574e-06, 0.4861839637949368, -0.021712454485320096, 9.470125832724226e-11, 8.959351751379364, 6.2831853071795845, 0.10000000178756306]
ramp_free_means_carbon=[0.2445537156141517, 0.07597303082211063,100.9519493198850647, 1.6017819143290766e-07, 3.315492244571699e-05, 0.08995451289980627, -0.003381307141896925, 1.625196327083214e-10, 8.959351751379364,  6.2831853071795845, 0.5415684133427566]
cmaes_ramped_time=ramp_fit.test_vals(ramp_means_carbon, likelihood="timeseries", test=False)
cmaes_rampfree_time=ramp_fit.test_vals(ramp_free_means_carbon, likelihood="timeseries", test=False)
cmaes_rampfree_time_2=ramp_fit.test_vals(ramp_free_means_black_2, likelihood="timeseries", test=False)

data_harmonics=harm_class.generate_harmonics(time_results, current_results)
ramp_harmonics=harm_class.generate_harmonics(time_results, cmaes_ramped_time)
ramp_free_harmonics=harm_class.generate_harmonics(time_results, cmaes_rampfree_time)
ramp_free_harmonics_2=harm_class.generate_harmonics(time_results, cmaes_rampfree_time_2)
#harm_class.plot_harmonics(time_results, method="abs", Experimental=data_harmonics, Ramped=ramp_harmonics)#, Ramp_free=ramp_free_harmonics)
harm_class.harmonics_and_time(ramp_fit.t_nondim(time_results), folder, "abs", \
                            Experimental_harmonics=ramp_fit.i_nondim(data_harmonics),Sinusoidal_harmonics=ramp_fit.i_nondim(ramp_free_harmonics),Ramped_harmonics=ramp_fit.i_nondim(ramp_harmonics),\
                             Experimental_time_series=ramp_fit.i_nondim(current_results),Sinusoidal_time_series=ramp_fit.i_nondim(cmaes_rampfree_time), Ramped_time_series=ramp_fit.i_nondim(cmaes_ramped_time), alpha=0.5) ##
test=ramp_fit.test_vals(ramp_means_black, likelihood="timeseries", test=False)
test_harmonics=harm_class.generate_harmonics(time_results, test)






param_boundaries=np.zeros((2, ramp_fit.n_parameters()))
for i in range(0, ramp_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[ramp_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[ramp_fit.optim_list[i]][1]

ramp_fit.define_boundaries(param_boundaries)

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
#found_parameters=x0
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
