import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_dcv import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]

folder_options=["/Black", "/Red", "/Plain", "/Carbon", "/Gold/Large", "/Gold/Small"]
dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/DcV/Black"#"/Black"#"/Gold/Large"#
Method ="PostScan"
type="asA_3"
path=dir_path+data_path+folder
files= os.listdir(path)
print files
for data in files:
    print data
    if (Method in data)  and (type in data):
        print Method
        file=open(path+"/"+data)
        results=np.loadtxt(file, skiprows=1)


time_results=results[:,0]
voltage_results=results[:,1]
current_results=results[:,2]
param_list={
    'E_start': min(voltage_results), #(starting dc voltage - V)
    'E_reverse': max(voltage_results),    #  (reverse dc voltage - V)
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
    'original_gamma': 1e-10,
    'k_0': 10000.0, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/2000),
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 50
}
plt.subplot(1,2,1)
plt.title("Current results")
plt.ylabel("Current")
plt.xlabel("time")
plt.plot(time_results, current_results)
plt.subplot(1,2,2)
plt.title("Input voltage")
plt.ylabel("voltage")
plt.xlabel("time")
plt.plot(time_results, voltage_results)
plt.show()
param_list['E_0']=0.23471918314326964
harmonic_range=np.arange(1,7,1)
dcv_fit=single_electron(param_list, params_for_opt, harmonic_range, 1.0)
dcv_fit.label="cmaes"
time_results=time_results[0::8]/dcv_fit.nd_param.c_T0
current_results=current_results[0::8]/dcv_fit.nd_param.c_I0
voltage_results=voltage_results[0::8]/dcv_fit.nd_param.c_E0
dcv_fit.time_vec=time_results
signal_length=len(current_results)
dcv_fit.num_points=signal_length
frequencies=np.fft.fftfreq(signal_length, dcv_fit.time_vec[1]-dcv_fit.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
dcv_fit.frequencies=frequencies
last_point= (harmonic_range[-1]*dcv_fit.nd_param.omega)+(dcv_fit.nd_param.omega*0.5)
plot_frequencies=frequencies[np.where(frequencies<last_point)]
dcv_fit.test_frequencies=plot_frequencies
harm_class=harmonics(harmonic_range, dcv_fit.nd_param.omega*dcv_fit.nd_param.c_T0, 0.05)
dcv_fit.label="MCMC"
dcv_fit.optim_list=[]#['E_0', 'k_0', 'Ru', 'Cdl','gamma', 'omega']
dcv_fit.pass_extra_data(current_results, False)
chains=np.load("ramped_results")
means=[np.mean(chains[:, 5000:, x]) for x in np.arange(0,len(dcv_fit.optim_list))]
dcv_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl','gamma', 'omega', 'phase', 'alpha']
param_bounds={
    'E_0':[0.1, 0.4],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 5e3],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-3], #(capacitance parameters)
    'CdlE1': [0,0.1],#0.000653657774506,
    'CdlE2': [0,0.1],#0.000245772700637,
    'CdlE3': [0,0.1],#1.10053945995e-06,
    'gamma': [1e-11,1e-9],
    'k_0': [0, 1e4], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    'phase' : [0, 2*math.pi]
}

dcv_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl',"CdlE1",'gamma', 'omega', 'phase', 'alpha']
param_boundaries=np.zeros((2, dcv_fit.n_parameters()))
for i in range(0, dcv_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[dcv_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[dcv_fit.optim_list[i]][1]

dcv_fit.define_boundaries(param_boundaries)
means=[1.92982653e-01, 3.52336088e+00, 4.99999999e+02, 4.46761269e-15, 1.30400052e-10, 8.89117171e+00]#830.276082338262
nums=1
rus=np.linspace(100, 5, nums)
means=[0.2501936847262157, 22.64822551278051, 1.4602089893583, 1e-4, 0.021593193177659842,5.8583709981749646e-11, 8.962676747582703, 6.2831853071777255, 0.5999999999991238]
gold_means=[0.2634471594256482, 182.46968143921333, 251.52483060070875, 9.948258528983337e-05, 0.021593193177659842, 6.486044409236416e-10, 8.941736782729983, 1.6346028726977888, 0.10000000000199746]
voltages=dcv_fit.define_voltages()
test=dcv_fit.simulate(means,frequencies, "yes", "timeseries", "no" )
def dim_current(current):
    return current*dcv_fit.nd_param.c_I0
def dim_time(time):
    return time*dcv_fit.nd_param.c_T0
def dim_potential(voltage):
    return voltage*dcv_fit.nd_param.c_E0

tr=((dcv_fit.nd_param.E_reverse-dcv_fit.nd_param.E_start));
dcv_fit.variable_returner()
means=[0.2501936847262157, 22.64822551278051, 1.4602089893583, 1e-4, 0.021593193177659842,5.8583709981749646e-11*0, 8.962676747582703, 6.2831853071777255, 0.5999999999991238]
test=dcv_fit.simulate(means,frequencies, "yes", "timeseries", "no" )
voltages=dim_potential(voltages)
#plt.plot(test)
plt.subplot(1,3,1)
plt.title("capactive current")
plt.plot(test)
plt.subplot(1,3,2)
means=[0.2501936847262157, 22.64822551278051, 1.4602089893583, 1e-4*0, 0.021593193177659842,5.8583709981749646e-11, 8.962676747582703, 6.2831853071777255, 0.5999999999991238]
test=dcv_fit.simulate(means,frequencies, "yes", "timeseries", "no" )
plt.title("faradaic current")
plt.plot(test)
plt.subplot(1,3,3)

plt.title("input potential")
plt.plot(voltages, test)

plt.show()

param_bounds={
    'E_0':[0.1, 0.4],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 3e3],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-2,2],#0.000653657774506,
    'CdlE2': [-0.1,0.1],#0.000245772700637,
    'CdlE3': [-0.05,0.05],#1.10053945995e-06,
    'gamma': [1e-11,4e-10],
    'k_0': [0, 1e4], #(reaction rate s-1)
    'alpha': [0.1, 0.9],
    'phase' : [0, 2*math.pi]
}
dcv_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl',"CdlE1", "CdlE2",'gamma', "alpha"]
param_boundaries=np.zeros((2, dcv_fit.n_parameters()))
for i in range(0, dcv_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[dcv_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[dcv_fit.optim_list[i]][1]
dcv_fit.define_boundaries(param_boundaries)
cmaes_problem=pints.SingleOutputProblem(dcv_fit, time_results, current_results)
score = pints.SumOfSquaresError(cmaes_problem)
dcv_fit.label="cmaes"
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(dcv_fit.optim_list))], [np.ones(len(dcv_fit.optim_list))])
x0=abs(np.random.rand(dcv_fit.n_parameters()))#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
print len(x0), dcv_fit.n_parameters()

found_parameters, found_value=pints.optimise(
                                            score,
                                            x0,
                                            boundaries=CMAES_boundaries,
                                            method=pints.CMAES
                                            )
cmaes_results=dcv_fit.change_norm_group(found_parameters, "un_norm")

test=dcv_fit.simulate(cmaes_results,frequencies, "no", "timeseries", "no" )
#dcv_fit.simulate(found_parameters,time_results, normalise=True, likelihood="fourier", test=True )
print cmaes_results
plt.plot(voltages, current_results)
plt.plot(voltages, test)
plt.show()
