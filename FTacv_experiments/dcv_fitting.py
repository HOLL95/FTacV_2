import numpy as np
import os
import matplotlib.pyplot as plt
#import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
#from single_e_class_dcv import single_electron
#from harmonics_plotter import harmonics
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

dec_amount=4
time_results=times[Method][0::dec_amount]
voltage_results=voltages[Method][0::dec_amount]
current_results=currents[Method][0::dec_amount]
plt.plot(voltage_results, current_results)
plt.show()
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
    'Cdlinv': 1e-6, #(capacitance parameters)
    'Cdlinv1': 0,#0.000653657774506,
    'Cdlinv2': 0,#0.000245772700637,
    'Cdlinv3': 0,#1.10053945995e-06,
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
plt.rcParams.update({'font.size': 20})
plt.subplot(1,2,2)
plt.title("DCV current vs voltage")
plt.ylabel("Current(mA)")
plt.xlabel("Voltage(V)")
plt.plot(voltage_results, current_results*1e3)
plt.subplot(1,2,1)
plt.title("DCV voltage vs time")
plt.ylabel("Voltage(V)")
plt.xlabel("Time(s)")
plt.plot(time_results, voltage_results)
plt.show()
param_list['E_0']=0.23471918314326964
harmonic_range=np.arange(1,7,1)
dcv_fit=single_electron(param_list, params_for_opt, harmonic_range, 1.0)
dcv_fit.label="cmaes"
tr=((dcv_fit.nd_param.E_reverse-dcv_fit.nd_param.E_start));
time_results=time_results/dcv_fit.nd_param.c_T0
current_results=current_results/dcv_fit.nd_param.c_I0
voltage_results=voltage_results/dcv_fit.nd_param.c_E0
time_idx=np.where(voltage_results==max(voltage_results))
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
    return current*dcv_fit.nd_param.c_I0
def dim_time(time):
    return time*dcv_fit.nd_param.c_T0
def dim_potential(voltage):
    return voltage*dcv_fit.nd_param.c_E0



param_bounds={
    'E_0':[0.1, 0.35],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [1, 1e3],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-1,1],#0.000653657774506,
    'CdlE2': [-0.1,0.1],#0.000245772700637,
    'CdlE3': [-0.03,0.03],#1.10053945995e-06,
    'Cdlinv': [0,1e-4], #(capacitance parameters)
    'Cdlinv1': [-1,1],#0.000653657774506,
    'Cdlinv2': [-0.1,0.1],#0.000245772700637,
    'Cdlinv3': [-0.03,0.03],#1.10053945995e-06,
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
dcv_fit.dispersion=True
dcv_fit.optim_list=["E0_mean", "E0_std","k_0","Ru", "Cdl","CdlE1","CdlE2","Cdlinv", "Cdlinv1","Cdlinv2", "gamma", "alpha"]
dcv_fit.dim_dict["CdlE3"]=0
dcv_fit.dim_dict["Cdlinv3"]=0
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
print(len(x0), dcv_fit.n_parameters())

found_parameters, found_value=pints.optimise(
                                            score,
                                            x0,
                                            boundaries=CMAES_boundaries,
                                            method=pints.CMAES
                                            )
cmaes_results=dcv_fit.change_norm_group(found_parameters, "un_norm")

test1=dcv_fit.simulate(cmaes_results,frequencies, "no", "timeseries", "no" )
#dcv_fit.simulate(found_parameters,time_results, normalise=True, likelihood="fourier", test=True )
print(list(cmaes_results))
plt.plot(dcv_fit.e_nondim(dcv_fit.voltages), dcv_fit.i_nondim(current_results), label="Data")
plt.plot(dcv_fit.e_nondim(dcv_fit.voltages), dcv_fit.i_nondim(test1), label="Simulation")
plt.legend()
plt.xlabel("Voltage(V)")
plt.ylabel("Current(A)")
plt.title(folder)
plt.show()
carbon_means=[ 3.22160701e-01,  2.40406450e+02,  8.97760719e-05,  3.29096335e-01,-5.10473887e-03,  2.05353151e-11]
