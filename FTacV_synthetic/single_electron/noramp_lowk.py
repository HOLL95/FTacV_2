import isolver
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class  import single_electron
import pints
import pints.plot
import copy
import time
from matplotlib.pyplot import figure
num_runs=1
filename="Convergence_results_"
points_range=int(1e4)
params_for_opt=['E_0', 'k_0', 'Ru', 'Cdl']
init_params=np.zeros((num_runs, len(params_for_opt)))
true_k=1e1
true_omega=32*math.pi
param_list={
    'E_start': -0.85, #(starting dc voltage - V)
    'E_reverse': -0.1,    #  (reverse dc voltage - V)
    'omega':true_omega,#8.88480830076,  #    (frequency Hz)
    'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 24.97e-3,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 20.0,  #     (uncompensated resistance ohms)
    'Cdl': 0.000134, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 6.5e-12,          # (surface coverage per unit area)
    'E_0': -0.4,      #       (reversible potential V)
    'E0_std':0.0312279186927,# (reversible potential dispersion)
    'k_0': true_k, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
    'time_end':1000,
    'num_peaks' : 30
}
harmonic_range=np.arange(1,13,1)
noramp=single_electron(param_list, params_for_opt, harmonic_range, 0.05)
freq_values=[20, 400, 60, 30, 500, 80, 1, 50, 70, 200, 100, 300, 40, 90, 10]
freq_values=np.sort(freq_values)
noramp.label="cmaes"
noramp.nd_param.time_end=(noramp.nd_param.num_peaks/noramp.nd_param.nd_omega)*2*math.pi#points_range*noramp.nd_param.sampling_freq
print noramp.nd_param.num_peaks/noramp.nd_param.omega
#noramp.nd_param.time_end=(2*(noramp.nd_param.E_reverse-noramp.nd_param.E_start)/noramp.nd_param.v)
noramp.times(points_range)
frequencies=np.fft.fftfreq(int(points_range), noramp.time_vec[1]-noramp.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
noramp.frequencies=frequencies
last_point= (harmonic_range[-1]*noramp.nd_param.omega)+(noramp.nd_param.omega*0.5)
plot_frequncies=frequencies[np.where(frequencies<last_point)]
noramp.test_frequencies=plot_frequncies
true_params=[-0.4, true_k, 6.5e-12,param_list['Cdl'], param_list['CdlE1'],param_list['CdlE2'],param_list['CdlE3']]
param_boundaries=[ [-0.6, 1, 1.5e-12,-1,0,-0.1,-0.1], \
                    [-0.2, 1e4, 1.5e-11,1,1,0.1,0.1]]# #
noramp.optim_list=[]
synthetic_data=noramp.simulate([], frequencies, "nah", "timeseries", "no")
plot_times=np.multiply(noramp.time_vec, noramp.nd_param.c_T0)
plot_current_data=np.multiply(synthetic_data, noramp.nd_param.c_I0)
figure(num=None, figsize=(8,5), dpi=120, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=0.15)
plt.plot(plot_times, plot_current_data)
plt.show()
#f=open('time_series_nr.fig', 'w')
#np.save(f, [plot_times, plot_current_data])
#f.close()
#plt.xlabel('Time(s)')
#plt.ylabel('Current(A)')
#plt.show()
noise_val=0.02
noise_max=max(synthetic_data)*noise_val
noise=np.random.normal(0,noise_max, len(synthetic_data))
noramp.define_boundaries(param_boundaries)
noisy_data=np.add(synthetic_data, noise)
noisy_data_filtered=noramp.kaiser_filter(noisy_data, False)
noramp.pass_extra_data(noisy_data, noisy_data_filtered)
noramp.simulate([], frequencies, "nah", "fourier", "yes")
noramp.optim_list=['E_0', 'k_0', 'Ru', 'Cdl']
cmaes_problem=pints.SingleOutputProblem(noramp, noramp.time_vec, noisy_data)
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
cmaes_runs=0
cmaes_results=np.zeros((cmaes_runs,len(param_boundaries[0])))
cmaes_scores=np.multiply(np.ones(cmaes_runs), 1e20);
k=0
found_parameters=np.zeros(len(param_boundaries[0]))
x0=abs(np.random.rand(len(param_boundaries[0])))
for i in range(0,len(true_params)):
    x0[i]=noramp.normalise(true_params[i], [param_boundaries[0][i],param_boundaries[1][i]])
while k<cmaes_runs:
    x0=abs(np.random.rand(len(param_boundaries[0])))
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
    cmaes_params=np.zeros(len(found_parameters))
    for i in range(0,len(found_parameters)):
        cmaes_params[i]=noramp.un_normalise(found_parameters[i], [param_boundaries[0][i],param_boundaries[1][i]])
    print cmaes_params
    cmaes_results[k,:]=cmaes_params
    cmaes_scores[k]=noramp.simulate(cmaes_params,frequencies, "Dont optimise", "fourier", "no", score=True )
    #noramp.pass_extra_data(noisy_data)
    #noramp.simulate(true_params, frequencies, "nah", "fourier", "yes")
    #noramp.pass_extra_data(noisy_data_filtered)
    k+=1


noramp.optim_list=['E_0', 'k_0', 'Cdl', 'Ru']
found_parameters=[param_list[x]for x in noramp.optim_list]
mcmc_problem=pints.SingleOutputProblem(noramp, np.linspace(0, 1, len(noisy_data_filtered)), noisy_data_filtered)
noramp.label="MCMC"
updated_lb=np.append(np.multiply([param_list[x]for x in noramp.optim_list], 0.5), [0*noise_max])#found_parameters[3]*0.97,
updated_ub=np.append(np.multiply([param_list[x]for x in noramp.optim_list], 1.5), [max(noisy_data_filtered)])#found_parameters[3]*1.03,
updated_boundaries=[updated_lb, updated_ub]
updated_boundaries=np.sort(updated_boundaries, 0)
noramp.define_boundaries(updated_boundaries)
print updated_boundaries
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                updated_boundaries[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
found_parameters=np.append(found_parameters, (70000))
print found_parameters
xs=[found_parameters,
    found_parameters,
    found_parameters
    ]
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()
rhats=pints.rhat_all_params(chains)
k_rhat=rhats[1]
#for k in range(0, len(chains)):
#    for i in range(0, len(noramp.optim_list)+1):
#        for j in range(0, len(chains[k,:,i])):
#            chains[k,j,i]=noramp.un_normalise(chains[k,j,i], [updated_lb[i],updated_ub[i]])
#k_means=np.zeros(3)
#e_means=np.zeros(3)
#for i in range(0, 3):
#    k_means[i]=np.mean(chains[i, 3000:, 1])
#    e_means[i]=np.mean(chains[i, 3000, 0])
#print k_means
#print np.mean(k_means)
#print k_rhat
#print "noise"+str(noise_val)
pints.plot.trace(chains[:, :, :])
print noise_max
plt.show()
save=raw_input("Do you want to save this data? Y/N?")
if save =="Y":
    filename=raw_input("filename?")
    f=open(filename, "w")
    np.save(f, chains)
    f.close()
#filename="populationMCMC_" + str(noise_val) + "_" + str(points_range[lcv_1])+ ".txt"
#f=open(filename, "w")
#np.save(f, chains)
#f.close()