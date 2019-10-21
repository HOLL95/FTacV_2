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
    'omega':8.977950e+00,#8.88480830076,  #    (frequency Hz)
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

param_list['E_0']=0.23471918314326964
harmonic_range=np.arange(4,7,1)
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

#param_boundaries=[[0.98*param_list['omega'], 0, 0.5], [0.98*param_list['omega'], 2*math.pi, 1.5]]
param_boundaries=[[0, 0.5], [2*math.pi, 1.5]]

ramp_fit.optim_list=[]
ramp_fit.pass_extra_data(ramp_fit.voltages, False)
test=ramp_fit.simulate([], frequencies, "no", "timeseries", "yes")
ramp_fit.pass_extra_data(test, False)
ramp_fit.optim_list=['phase', 'v']
ramp_fit.define_boundaries(param_boundaries)
cmaes_problem=pints.SingleOutputProblem(ramp_fit, time_results, test)
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
x0=np.ones(len(ramp_fit.optim_list))*0.15
found_parameters, found_value=pints.optimise(
                                            score,
                                            x0,
                                            boundaries=CMAES_boundaries,
                                            method=pints.CMAES
                                            )
ramp_fit.simulate(found_parameters, frequencies, "optimise", "timeseries", "yes")
print found_parameters
