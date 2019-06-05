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
harmonic_range=np.arange(2,9,1)
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
ramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma', 'omega', 'phase','alpha']

means=[1.92982653e-01, 3.52336088e+00, 4.99999999e+02, 4.46761269e-15, 1.30400052e-10, 8.89117171e+00]#830.276082338262
nums=1
rus=np.linspace(100, 5, nums)

#means=[2.41444794e-01, 6.17012486e-01, 1.04871338e+01, 3.51866140e-06, 4.49820672e-10, 8.94073199e+00, 1.39984093e-01]#, 4.56684293e+00]
means=[0.23625846615525375, 282.2528765716474, 150.32008743633, 7.233724535557822e-06, 6.1720248951887567e-11, 8.940733606629868, 1.4114762214007537, 0.4]
#means=[0.2724040125447893, 121.81748721519133, 1207.1166997969588, 5.295143515633727e-06, 1.1895825195613567e-10, 8.940848161584075, 1.8671427717422973, 0.10000000005872724]
test=ramp_fit.simulate(means,frequencies, "no", "timeseries", "no" )

plt.show()

plt.plot(current_results)
plt.plot(np.array(test)*-1)
plt.show()
noise_p=0.02
noise_level=max(test)*noise_p
noise=np.random.normal(0, noise_level, len(test))
synthetic_data=np.add(noise, test)
fourier_test=ramp_fit.kaiser_filter(synthetic_data)
test_harmonics=harm_class.generate_harmonics(time_results, test)
data_harmonics=harm_class.generate_harmonics(time_results, current_results)
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)
plt.plot(time_results, test, label="numerical")
plt.plot(time_results, current_results, alpha=0.7)
plt.legend()
plt.show()
harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics,"abs", "numerical", "data")
likelihood_func=ramp_fit.kaiser_filter(current_results)
param_bounds={
    'E_0':[0.1, 0.4],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 1e4],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-3], #(capacitance parameters)
    'CdlE1': [0,0.1],#0.000653657774506,
    'CdlE2': [0,0.1],#0.000245772700637,
    'CdlE3': [0,0.1],#1.10053945995e-06,
    'gamma': [1e-11,1e-9],
    'k_0': [0, 1e4], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    'phase' : [0, 2*math.pi]
}

ramp_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl',"CdlE1",'gamma', 'omega', 'phase', 'alpha']
param_boundaries=np.zeros((2, ramp_fit.n_parameters()))
for i in range(0, ramp_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[ramp_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[ramp_fit.optim_list[i]][1]

dummy_times=np.linspace(0, 1, len(likelihood_func))
ramp_fit.define_boundaries(param_boundaries)
cmaes_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
ramp_fit.label="cmaes"
#ramp_fit.method_label="timeseries"
#cmaes_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, fourier_test)
plt.plot(test)
plt.show()
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
x0=abs(np.random.rand(ramp_fit.n_parameters()))#ramp_fit.change_norm_group(means[1:], "norm")#abs(np.random.rand(ramp_fit.n_parameters()))
ramp_fit.pass_extra_data(current_results, likelihood_func)
for i in range(0, 1):
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
    x0=found_parameters
    print found_parameters

ramp_fit.simulate(found_parameters, frequencies, "optimise", "fourier", "yes")
ramp_fit.simulate(found_parameters, frequencies, "optimise", "timeseries", "yes")
cmaes_results=ramp_fit.change_norm_group(found_parameters, "un_norm")
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
