
import isolver
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class  import single_electron
import pints
import pints.plot
import copy
import time
num_runs=1
filename="Convergence_results_"
points_range=int(1e6)
params_for_opt=['E_0', 'k_0', 'Ru', 'Cdl']
init_params=np.zeros((num_runs, len(params_for_opt)))
true_k=1e4
true_omega=8.8
param_list={
    'E_start': -0.85, #(starting dc voltage - V)
    'E_reverse': -0.1,    #  (reverse dc voltage - V)
    'omega':true_omega,#8.88480830076,  #    (frequency Hz)
    'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 24.97e-3,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 20.0,  #     (uncompensated resistance ohms)
    'Cdl':3.7e-5, #(capacitance parameters)
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
    'time_end':1000
}
harmonic_range=np.arange(4,9,1)
noramp=single_electron(param_list, params_for_opt, harmonic_range, 0.5)
noramp.label="cmaes"
#noramp.nd_param.time_end=points_range*noramp.nd_param.sampling_freq
noramp.nd_param.time_end=2*((noramp.nd_param.E_reverse-noramp.nd_param.E_start)/noramp.nd_param.v)
noramp.times(points_range)
#noramp.times(int(noramp.nd_param.time_end/noramp.nd_param.sampling_freq))
#points_range=len(noramp.time_vec)
frequencies=np.fft.fftfreq(int(points_range), noramp.time_vec[1]-noramp.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
noramp.frequencies=frequencies
#(frequencies>((harmonic_range[0]*noramp.nd_param.omega)-(noramp.nd_param.omega*0.5))))
frequencies=frequencies[np.where((frequencies>((harmonic_range[0]*noramp.nd_param.omega)-(noramp.nd_param.omega*0.5)))&(frequencies<((harmonic_range[-1]*noramp.nd_param.omega)+(noramp.nd_param.omega*0.5))))]
noramp.test_frequencies=frequencies

true_params=[-0.4, true_k, 6.5e-12,param_list['Cdl'], param_list['CdlE1'],param_list['CdlE2'],param_list['CdlE3']]
param_boundaries=[ [-0.6, 1, 1.5e-12,-1,0,-0.1,-0.1], \
                    [-0.2, 1e4, 1.5e-11,1,1,0.1,0.1]]# #
noramp.optim_list=[]
synthetic_data=noramp.simulate([], frequencies, "nah", "timeseries")
noise_val=0.0
noise_max=noise_val*max(synthetic_data)
noise=np.random.normal(0,noise_max, len(synthetic_data))
noramp.define_boundaries(param_boundaries)
Cdl_range=[2*math.pi, 16*math.pi, 64*math.pi]
#noramp.optim_list=['E_0', 'k_0', 'Ru', 'Cdl']
noramp.optim_list=[ 'Cdl']
f=np.fft.fftfreq(int(points_range), noramp.time_vec[1]-noramp.time_vec[0])
L=int(points_range)
window=np.hanning(L)
pi=math.pi
Cdl_range=[1e-5, 1e-4, 3.7e-3]#[4*pi, 16*pi,  64*pi, 128*pi, 256*math.pi]
fig, ax=plt.subplots(noramp.num_harmonics,1)
for i in range(0,len(Cdl_range)):
    synthetic_data=noramp.simulate([Cdl_range[i]], frequencies, "nah", "timeseries")
    noisy_data=np.add(synthetic_data,noise)
    harmonics=noramp.kaiser_filter(noisy_data, True)
    print len(harmonics)
    for j in range(0, noramp.num_harmonics):
        ax[j].plot(harmonics[j,:], label=noramp.optim_list[0]+" " + str(Cdl_range[i]), alpha=0.5)
    #plt.show()
    #Y=np.fft.fft(np.multiply(noisy_data, window))
    #if i==len(Cdl_range):
    #    plt.plot(f[0:len(f)/2],Y[0:len(Y)/2], label=("Cdl "+str(Cdl_range[i])), alpha=0.3)
    #else:
    #    plt.plot(f[0:len(f)/2],Y[0:len(Y)/2], label=("Cdl "+str(Cdl_range[i])))
    #plt.legend()
#plt.show()
#plt.plot(synthetic_data)
#plt.show()
plt.legend()
plt.show()

#noisy_data_filtered=noramp.kaiser_filter(noisy_data, frequencies)

#noramp.pass_extra_data(noisy_data, noisy_data_filtered)
"""
noramp.simulate(true_params, frequencies, "nah", "fourier", "no")
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


#z=np.where(cmaes_scores==min(cmaes_scores))
#found_parameters=cmaes_results[z[0],:][0]

noramp.optim_list=['E_0', 'k_0', 'Cdl', 'Ru']
found_parameters=[param_list[x]for x in noramp.optim_list]
mcmc_problem=pints.SingleOutputProblem(noramp, noramp.time_vec, noisy_data)
noramp.label="MCMC"
updated_lb=np.append(np.multiply([param_list[x]for x in noramp.optim_list], 0.8), [0*noise_max])#found_parameters[3]*0.97,
updated_ub=np.append(np.multiply([param_list[x]for x in noramp.optim_list], 1.2), [2*noise_max])#found_parameters[3]*1.03,
updated_boundaries=[updated_lb, updated_ub]
updated_boundaries=np.sort(updated_boundaries, 0)
noramp.define_boundaries(updated_boundaries)
print updated_boundaries
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                updated_boundaries[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
found_parameters=np.append(found_parameters, noise_max)
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
"""
