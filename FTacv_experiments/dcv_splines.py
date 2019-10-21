import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import scipy
import matplotlib.pyplot as plt
from single_e_class_dcv import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]
def file_opener(files, file_path, file_dict):
    voltage_results={}
    time_results={}
    current_results={}
    experiments=list(file_dict.keys())
    for data in files:
        for keys in experiments:
            Method=keys
            type=file_dict[keys]
            if (Method in data)  and (type in data):
                file=open(path+"/"+data)
                results=np.loadtxt(file, skiprows=1)
                voltage_results[Method]=results[:,1]
                time_results[Method]=results[:,0]
                current_results[Method]=results[:,2]
                break
    return voltage_results, current_results, time_results
folder_options=["/Black", "/Red", "/Plain", "/Carbon", "/Gold/Large", "/Gold/Small"]
dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/DcV/Carbon"#"/Black"#"/Gold/Large"#
Method ="PostSinusoidal"
type="asA_4"
path=dir_path+data_path+folder
files= os.listdir(path)
file_dict={"PostSinusoidal":"asA_4", "IntScan":"asA_4", "PostLinearRamp":"asA_4"}
voltages, currents, times=file_opener(files, path, file_dict)
for files in list(file_dict.keys()):
    plt.plot(voltages[files], currents[files], label=files)
plt.legend()
plt.show()


time_results=times[Method]
voltage_results=voltages[Method]
current_results=currents[Method]

param_list={
    "E_0":0.2,
    'E_start': min(voltage_results), #(starting dc voltage - V)
    'E_reverse': max(voltage_results),    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': max(voltage_results)-min(voltage_results),   #(ac voltage amplitude - V) freq_range[j],#
    'v': 29.8e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 20.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-6, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,          # (surface coverage per unit area)
    'original_gamma': 1e-10,
    "cap_phase":0,
    'k_0': 10000.0, #(reaction rate s-1)
    "E0_mean":0.2,
    "E0_std": 0.09,
    "k0_shape":0.954,
    "k0_loc":100,
    "k0_scale":50,
    "k0_range":1e3,
    "cap_phase":0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 50
}

harmonic_range=np.arange(1,7,1)
dcv_fit=single_electron(param_list, params_for_opt, harmonic_range, 1.0)
dcv_fit.label="cmaes"
tr=((dcv_fit.nd_param.E_reverse-dcv_fit.nd_param.E_start));
time_results=time_results/dcv_fit.nd_param.c_T0
current_results=current_results/dcv_fit.nd_param.c_I0
voltage_results=voltage_results/dcv_fit.nd_param.c_E0
#time_idx=np.where(voltage_results==max(voltage_results))
#time_results=time_results[:time_idx[0][0]]
#voltage_results=voltage_results[:time_idx[0][0]]
#current_results=current_results[:time_idx[0][0]]
dcv_fit.voltages=voltage_results
dcv_fit.initial_val=current_results[0]
print(time_results[1]-time_results[0])
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
dcv_fit.dispersion=False
means=[param_list[dcv_fit.optim_list[x]] for x in np.arange(0,len(dcv_fit.optim_list))]
def dim_current(current):
    return np.multiply(current, dcv_fit.nd_param.c_I0)
def dim_time(time):
    return np.multiply(time, dcv_fit.nd_param.c_T0)
def dim_potential(voltage):
    return np.multiply(voltage, dcv_fit.nd_param.c_E0)
param_bounds={
    'E_0':[0.12, 0.4],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [1, 1e3],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-1,1],#0.000653657774506,
    'CdlE2': [-0.1,0.1],#0.000245772700637,
    'CdlE3': [-0.005,0.005],#1.10053945995e-06,
    'gamma': [1e-11,4e-10],
    'k_0': [1, 1e2], #(reaction rate s-1)
    'alpha': [0.1, 0.9],
    'phase' : [0, 2*math.pi],
    "E0_mean":[0.2, 0.4],
    "E0_std": [0, 0.2],
    "k0_shape":[0,5],
    "k0_loc":[1, 1e4],
    "k0_scale":[0,2e3],
    "k0_range":[1e2, 1e4],
}
scan_params=['E_0','k_0', "Ru", "Cdl","CdlE1","CdlE2", "gamma"]
carbon_means=[0.2712720627064147, 53.00972798907347, 126.43581153952566, 3.160962997490096e-05, 0.14647947194537103, -0.005848319334033306, 1.0072445202882476e-10, 0]
dcv_fit.dispersion=True
scan_params=[]
for q in range(0, len(scan_params)):
    dcv_fit.dim_dict[scan_params[q]]=carbon_means[q]
for i in range(0, len(scan_params)):
    plt.subplot(2,4,i+1)
    plt.title(scan_params[i])
    param_range=np.linspace(param_bounds[scan_params[i]][0], param_bounds[scan_params[i]][1], 4)
    dcv_fit.optim_list=[scan_params[i]]
    for j in range(0, 4):
        print(dcv_fit.nd_param.CdlE1)
        test=dcv_fit.simulate([param_range[j]],frequencies, "yes", "timeseries", "no" )
        plt.plot(dcv_fit.e_nondim(voltage_results), dcv_fit.i_nondim(test), label=param_range[j])
    plt.plot(dcv_fit.e_nondim(voltage_results),dcv_fit.i_nondim(current_results), color="black")
    plt.legend()
    for q in range(0, len(scan_params)):
        dcv_fit.dim_dict[scan_params[q]]=carbon_means[q]
plt.show()

def nan_helper(y):
 return np.isnan(y), lambda z: z.nonzero()[0]
dcv_fit.optim_list=["E0_mean", "E0_std", "k_0", "Ru", "alpha"]
test=dcv_fit.simulate([0.23264569588050996, 0.0674603084243984, 84.4445611107851, 227.55356998656478, 0.5], frequencies)
test2=dcv_fit.simulate([0.26309322667095436, 0.08411919477796326, 999.9999999749896, 5.005163812028062e-10, 0.7999999999895195], frequencies)
test3=dcv_fit.simulate([0.24042307692307694, 0.05995192307692309, 81.25609725526283, 100.0,0.6], frequencies)
test4=dcv_fit.simulate([0.2507124585192858, 0.012810161448688611, 81.25609725526283,100, 0.7811192097028208], frequencies)
subtract=np.subtract(current_results, test)
#plt.plot(voltage_results, subtract)
#plt.plot(voltage_results, current_results)
#plt.show()
f_voltage_range=np.divide(param_bounds["E_0"], dcv_fit.nd_param.c_E0)
f_voltage_idx=np.where((voltage_results>f_voltage_range[0]) & (voltage_results<f_voltage_range[1]))
cap_voltages=voltage_results[f_voltage_idx]
cap_current=copy.deepcopy(current_results)#np.subtract(current_results, test)
cap_current[f_voltage_idx]=float('NaN')
nans, x=nan_helper(cap_current)
cap_current[nans]=np.interp(x(nans), x(~nans), cap_current[~nans])
max_idx=np.where(voltage_results==(max(voltage_results)))
max_idx=max_idx[0][0]
print(max_idx)
section_1=cap_current[:max_idx]
section_2=cap_current[max_idx:]
degree=3

section_1_polyfit=np.poly1d(np.polyfit(time_results[:max_idx], section_1, degree))
section_2_polyfit=np.poly1d(np.polyfit(time_results[max_idx:], section_2,degree))
total_polyfit= np.append(section_1_polyfit(time_results[:max_idx]), section_2_polyfit(time_results[max_idx:]))
cap_current_2=copy.deepcopy(cap_current)
cap_current[nans]=total_polyfit[nans]
cap_time=copy.deepcopy(time_results)
plt.plot(voltage_results, current_results)
plt.plot(voltage_results, np.add(cap_current, test2))
plt.plot(voltage_results, np.add(cap_current, test3))
plt.plot(voltage_results, np.add(cap_current, test4))
plt.show()


coarse_time=time_results[0::64]
coarse_current=current_results[0::64]
coarse_voltage=voltage_results[0::64]
coarse_time=np.linspace(0, time_results[-1], 18)
coarse_current=np.interp(coarse_time, time_results, current_results)
coarse_voltage=np.interp(coarse_time, time_results, voltage_results)
cs=scipy.interpolate.CubicSpline(coarse_time, coarse_current)
cs_results=cs(coarse_time)
plt.plot(coarse_voltage, coarse_current, "x")
plt.plot(coarse_voltage, cs_results)
#plt.plot(voltage_results, np.add(cs_results, test))
#plt.plot(voltage_results, np.add(cs_results, test2))
plt.show()


param_boundaries=np.zeros((2, dcv_fit.n_parameters()))
for i in range(0, dcv_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[dcv_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[dcv_fit.optim_list[i]][1]
dcv_fit.define_boundaries(param_boundaries)
dcv_fit.label="cmaes"
x0=abs(np.random.rand(dcv_fit.n_parameters()))#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
print(len(x0), dcv_fit.n_parameters())
