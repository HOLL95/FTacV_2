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
    'Ru': 300.0 ,  #     (uncompensated resistance ohms)
    'Cdl': 0.00000134, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1.4069958184991404e-10,          # (surface coverage per unit area)
    'k_0': 3.33567800e+03, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/100),
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 50
    }

param_list['E_0']=0.23471918314326964
harmonic_range=np.arange(4,7,1)
ramp_fit=single_electron(param_list, params_for_opt, harmonic_range, 1.0)
ramp_fit.label="cmaes"
ramp_fit.voltages=voltages[0::dec_amount, 1]/ramp_fit.nd_param.c_E0
time_results=time_results/ramp_fit.nd_param.c_T0
print time_results
current_results=current_results/ramp_fit.nd_param.c_I0
#plt.plot(time_results, current_results, label="data")
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
ramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma', 'omega']
ramp_fit.pass_extra_data(current_results, False)
chains=np.load("ramped_results")
means=[np.mean(chains[:, 5000:, x]) for x in np.arange(0,len(ramp_fit.optim_list))]
means=[2.3264595e-01, 2.51068461e+01, 3.80457587e+02, 7.27571754e-04,1.32687462e-10, 8.94056135e+00]

data_harmonics=harm_class.generate_harmonics(time_results, current_results)
for i in range(0, 1):#len(param_guesses)):
    time_series=ramp_fit.simulate(means, time_results, "no", "timeseries", "no")
    plt.plot(time_results, time_series)
    plt.plot(time_results, current_results, alpha=0.7)
    plt.show()
    exp_harmonics=harm_class.generate_harmonics(time_results, time_series)
    harm_class.plot_harmonics(time_results, exp_harmonics, data_harmonics)





"""
param_guesses=np.array([
[2.40666629e-01, 1.33659203e+04, 8.32138327e+00, 6.48986896e-0,  6.99137283e-10 ],
 [2.89424312e-01, 1.42829102e+04, 1.16728294e+00, 2.95794563e-02, 1.54470074e-10 ],
 [2.37426876e-01, 1.40219177e+04, 1.19240444e+01, 2.29657002e-01, 1.42238592e-09 ],
 [2.41984331e-01, 1.03875610e+00, 2.04067494e+01, 5.02952291e-03, 3.85177644e-10 ],
[2.34130799e-01, 4.02039142e+01, 2.46522095e+02, 3.91211938e-04, 2.59101097e-10 ],
 [2.77009135e-01, 1.01677918e+04, 1.22365956e+01, 1.17780590e+00, 1.21045676e-09 ],
 [2.85687934e-01, 1.39493817e+04, 5.61819373e-01, 4.73712354e-02, 1.46875670e-10 ],
 [2.85658802e-01, 1.31128704e+04, 5.60503182e-01, 4.73852493e-02, 1.46876400e-10 ],
 [2.74746341e-01, 1.57993170e+03, 1.46569159e+02, 3.42081253e-03, 6.72772614e-10],
 [2.84937019e-01, 1.24812926e+04, 1.51167967e+02, 1.58554929e-04, 1.47501166e-10]
])
#plt.show()
"""
