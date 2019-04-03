import sys
import time
import matplotlib.pyplot as plt
from params_class import params
from multiple_e_class import multiple_electron
import math
import numpy as np
import pints
import pints.plot
num_species=4
param_dict={
    'E_start': 0.6, #(starting dc voltage - V)
    'E_reverse':-0.1,    #  (reverse dc voltage - V)
    'omega':8.8,  #    (frequency Hz)
    'd_E': 20e-3,  #(ac voltage amplitude - V) freq_range[j],#
    'v': -0.1043081,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 20.0,  #     (uncompensated resistance ohms)
    'Cdl': 3.7e-3,#0.000133878548046, #(capacitance parameters)
    'CdlE1':  0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 6.5e-12,          # (surface coverage per unit area)
    'E_0': np.linspace(0.5, 0.2, num_species),   #       (reversible potential V),-0.016]
    'E0_std':0.0312279186927,# (reversible potential dispersion)
    'k_0': abs(np.multiply(np.random.rand(num_species), 1e1)), #(reaction rate s-1)
    'num_species': num_species,
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
}
method="classical"
seq_e=multiple_electron(param_dict,method , np.arange(3, 9))
if method=="noramp":
    time_end=2e4*seq_e.nd_param.sampling_freq
elif method == "classical":
    time_end=abs(2*((seq_e.nd_param.E_reverse-seq_e.nd_param.E_start)/seq_e.nd_param.v))
params_for_opt=['E_0', 'k_0']#, 'gamma', 'omega','Cdl','CdlE1','CdlE2','CdlE3']
seq_e.optimisation_params(params_for_opt)
times=seq_e.time_params(time_end)
l_frequencies=np.fft.fftfreq(len(times), param_dict['sampling_freq'])
end=np.where(l_frequencies>seq_e.harmonic_range[-1])
frequencies=l_frequencies[0:end[0][0]]
true_data=seq_e.simulate([], frequencies, gen_data=True)
plt.plot(true_data)
plt.show()
window=np.hanning(len(true_data))
Y=np.fft.fft(np.multiply(true_data, window))
plt.plot(l_frequencies, Y)
plt.show()
noise_val=0.02
noise_max=noise_val*max(true_data)
noise=np.random.normal(1,noise_max, len(true_data))
noise_mean=np.mean(noise)
noise_std=np.std(noise)
noisy_data=np.add(true_data,noise)
plt.plot(l_frequencies, np.fft.fft(noisy_data))
plt.show()
data_filtered=seq_e.top_hat_filter(frequencies, true_data)
noisy_data_filtered=seq_e.top_hat_filter(frequencies, noisy_data)
plt.plot(frequencies, data_filtered)
plt.plot(frequencies, noisy_data_filtered, alpha=0.7)
plt.show()
seq_e.pass_extra_data(noisy_data_filtered)
e0_lower_bound=[param_dict['E_0'][-1]]*num_species
e0_upper_bound=[param_dict['E_0'][0]]*num_species
k0_lower_bound=[10]*num_species
k0_upper_bound=[1e2]*num_species
k0_lower_bound[0]=param_dict['k_0'][0]
k0_upper_bound[0]=param_dict['k_0'][0]
other_params_lower_bound=[1.5e-12, 0,-1,-0.1,-0.1]
other_params_upper_bound=[1.5e-11, 1,1,0.1,0.1]
lower_bound=np.concatenate([e0_lower_bound, k0_lower_bound, [0]])#, other_params_lower_bound])
upper_bound=np.concatenate([e0_upper_bound, k0_upper_bound, [100]])#, other_params_upper_bound])
boundaries=[lower_bound, upper_bound]
seq_e.define_boundaries(boundaries)
seq_e.optim_list=['E_0', 'k_0']
print len(boundaries[0])
print [np.zeros(len(boundaries[0]))], [np.ones(len(boundaries[0]))]
"""
cmaes_problem=pints.SingleOutputProblem(seq_e, frequencies, noisy_data_filtered)
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.Boundaries([np.zeros(len(boundaries[0]))], [np.ones(len(boundaries[0]))])
for i in range(0,5):
    x0=abs(np.random.rand(len(boundaries[0])))
    print x0
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
    cmaes_params=np.zeros(len(found_parameters))
    for i in range(0,len(found_parameters)):
        cmaes_params[i]=seq_e.un_normalise(found_parameters[i], [boundaries[0][i],boundaries[1][i]])
    for i in range(0, len(params_for_opt)):
        print param_dict[params_for_opt[i]]
    print cmaes_params
    #seq_e.simulate(cmaes_params, frequencies, gen_data=False, test=True)
"""
mcmc_problem=pints.SingleOutputProblem(seq_e, frequencies, noisy_data_filtered)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(np.zeros(len(boundaries[0])),
                                np.ones(len(boundaries[0])))
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
initial_pos=np.ones(len(boundaries[0]))
xs=[initial_pos*0.5,
    initial_pos*0.5,
    initial_pos*0.5,
    ]
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.PopulationMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()

for k in range(0, len(chains)):
    for i in range(0, len(seq_e.optim_list)):
        for j in range(0, len(chains[k,:,i])):
            chains[k,j,i]=seq_e.un_normalise(chains[k,j,i], [boundaries[0][i],boundaries[1][i]])
pints.plot.trace(chains[:, 1000:, :])
plt.show()
