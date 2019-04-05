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
Method ="O_"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        voltages=np.loadtxt(path+"/"+data)
dec_amount=32
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
"""
param_list={
    'E_start': -180e-3, #(starting dc voltage - V)
    'E_reverse': 620e-3,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 29.8e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 300.0 ,  #     (uncompensated resistance ohms)
    'Cdl': 0.00534, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,          # (surface coverage per unit area)
    'k_0': 3.33567800e+03, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 50
    }
"""
de=300e-3
estart=260e-3-de
ereverse=estart+2*de

param_list={
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
    'gamma': 1e-10,          # (surface coverage per unit area)
    'k_0': 10000.0, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/2000),
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 50
}

param_list['E_0']=0.23471918314326964
harmonic_range=np.arange(3,7,1)
ramp_fit=single_electron(param_list, params_for_opt, harmonic_range, 1.0)
ramp_fit.label="cmaes"
ramp_fit.voltages=voltages[0::dec_amount, 1]/ramp_fit.nd_param.c_E0
time_results=time_results/ramp_fit.nd_param.c_T0
print time_results
current_results=current_results/ramp_fit.nd_param.c_I0
#plt.plot(time_results, ramp_fit.voltages, label="data")
ramp_fit.time_vec=time_results
signal_length=len(current_results)
ramp_fit.num_points=signal_length
frequencies=np.fft.fftfreq(signal_length, ramp_fit.time_vec[1]-ramp_fit.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
ramp_fit.frequencies=frequencies
last_point= (harmonic_range[-1]*ramp_fit.nd_param.omega)+(ramp_fit.nd_param.omega*0.5)
plot_frequencies=frequencies[np.where(frequencies<last_point)]
ramp_fit.test_frequencies=plot_frequencies
harm_class=harmonics(harmonic_range, ramp_fit.nd_param.omega*ramp_fit.nd_param.c_T0, 0.05)
ramp_fit.label="MCMC"
ramp_fit.optim_list=[]#['E_0', 'k_0', 'Ru', 'Cdl','gamma', 'omega']
ramp_fit.pass_extra_data(current_results, False)
chains=np.load("ramped_results")
means=[np.mean(chains[:, 5000:, x]) for x in np.arange(0,len(ramp_fit.optim_list))]
ramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma', 'omega']#, 'phase']
means=[0.233, 120.57601349e+00, 2.83225263e+02, 8.70271625e-04, 3.78984258e-10, 8.94070468e+00]
#means=[4.53228947e-01, 1.46151169e+02, 5.92813367e+02, 1.21792423e-04,3.10185626e-10, 8.94071394e+00]
#means=[2.86965303e-01, 1.21262967e+01, 3.62505000e+02, 1.29491311e-04,1.35050241e-10, 8.94070555e+00]
means=[2.33152423e-01, 45.83534674e+00, 2.04498067e+02, 7.53421861e-04, 1.17532353e-10, 8.94050130e+00]#, 8.07572569e-01]
#means=[2.85982586e-01, 5.57540878e+00, 3.62827219e+00, 1.86581548e-06, 1.20758344e-10, 8.94056256e+00, 3.64861364e+00]

#means=[2.86821677e-01, 1.23819142e+02, 3.65239741e+02, 3.30105547e-04, 1.34921412e-10, 8.94069398e+00]

test=ramp_fit.simulate(means,frequencies, "no", "timeseries", "no" )
plt.plot(test)
plt.show()
test_harmonics=harm_class.generate_harmonics(time_results, test)
data_harmonics=harm_class.generate_harmonics(time_results, current_results)
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)
plt.plot(time_results, test, label="numerical")
plt.legend()
plt.show()
harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics,"abs", "numerical", "data")
likelihood_func=ramp_fit.kaiser_filter(current_results)
param_boundaries=[[param_list['E_start'],1, 0, 0, 5.0e-11, 0.99*param_list['omega']], \
                    [param_list['E_reverse'], 1e4,500, 0.1, 5e-10, 1.01*param_list['omega']]]# #
ramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma', 'omega']
dummy_times=np.linspace(0, 1, len(likelihood_func))
ramp_fit.define_boundaries(param_boundaries)
#cmaes_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
ramp_fit.label="cmaes"
#ramp_fit.method_label="timeseries"
cmaes_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, likelihood_func)
plt.plot(test)
plt.show()
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
x0=ramp_fit.change_norm_group(means, "norm")

for i in range(0, 1):
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
    x0=found_parameters
    print found_parameters
cmaes_results=ramp_fit.change_norm_group(found_parameters, "un_norm")
print cmaes_results
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
