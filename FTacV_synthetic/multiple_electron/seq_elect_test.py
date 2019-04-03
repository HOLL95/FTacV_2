import sys
import time
import matplotlib.pyplot as plt
from params_class import params
import math
import numpy as np

num_species=3
param_list={),
    'E_start': 0.6, #(starting dc voltage - V)
    'E_reverse':-0.1,    #  (reverse dc voltage - V)
    'omega':60.05,  #    (frequency Hz)
    'd_E': 20e-3,  #(ac voltage amplitude - V) freq_range[j],#
    'v': -0.1043081,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 20.0,  #     (uncompensated resistance ohms)
    'Cdl': 0.000008,#0.000133878548046, #(capacitance parameters)
    'CdlE1':  0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 6.5e-12,          # (surface coverage per unit area)
    'E_0': np.linspace(0.6, -0.1, num_species),   #       (reversible potential V),-0.016]
    'E0_std':0.0312279186927,# (reversible potential dispersion)
    'k_0': abs(np.multiply(np.random.rand(num_species), 1e3)), #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 0#3*(math.pi/2),
}
key_list=param_list.keys()
nd_param=params(param_list)
for i in range(0, len(key_list)):
    nd_param.non_dimensionalise(key_list[i], param_list[key_list[i]])
if sys.argv[1]=="noramp":
    import seq_electron_noramp as seq_electron
    duration = int(1e4)
    nd_param.phase=3*(math.pi/2)
elif sys.argv[1]=="classical":
    import seq_electron_classical as seq_electron
    duration=int((2*((nd_param.E_reverse-nd_param.E_start)/nd_param.v))/nd_param.sampling_freq)
num_species=len(nd_param.E_0)
start= time.time()
print nd_param.E_0
print nd_param.E_start, nd_param.E_reverse
results=seq_electron.current_solver(nd_param.Cdl,nd_param.CdlE1,nd_param.CdlE2,nd_param.CdlE3,nd_param.nd_omega,nd_param.v, nd_param.alpha ,
                                    nd_param.E_start,  nd_param.E_reverse,  nd_param.d_E,  nd_param.Ru, nd_param.sampling_freq, duration,  nd_param.gamma, \
                                     nd_param.phase,num_species, nd_param.k_0, nd_param.E_0)
print time.time()-start
plt.plot(results[1][10:], results[0][10:])
plt.show()
window=np.hanning(len(results[1]))
time_series=np.multiply(window, results[0])
f=np.fft.fftfreq(len(results[1]), param_list['sampling_freq'])
plt.plot(f,np.fft.fft(time_series))
plt.show()
