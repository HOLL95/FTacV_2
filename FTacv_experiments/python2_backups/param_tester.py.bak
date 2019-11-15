import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_noramp  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]
dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Black"#"/Gold/Large"#
Method ="N_"#"GoldLarge_1"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        results2=np.loadtxt(path+"/"+data)
dec_amount=8
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=results2[0::dec_amount, 1]

current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=results2[0::dec_amount, 1]

length_list=[1e4, 2e4, 3e4]
dec_list=[8, 16, 32, 64]
repeat_num=5
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.1, #(electrode surface area cm^2)
    'Ru': 382,  #     (uncompensated resistance ohms)
    'Cdl': 5.440677575193328e-05, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1.9e-10,          # (surface coverage per unit area)
    'k_0': 357, #(reaction rate s-1)
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 1.406746302896052,
    "time_end": None,
    'num_peaks':2
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":False,
    "test": False,
    "likelihood":likelihood_options[0],
    "numerical_method": solver_list[3],
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,9,1),
    "experiment_time": None,
    "experiment_current": None,
    "experiment_voltage":None,
    "bounds_val":2000,
    "signal_length":int(2e4),
}
param_list['E_0']=0.217740555023939#(param_list['E_reverse']-param_list['E_start'])/2
noramp_params=single_electron(param_list, simulation_options, other_values)
harm_class=harmonics(other_values["harmonic_range"], noramp_params.nd_param.omega*noramp_params.nd_param.c_T0, 0.5)
noramp_params.label="MCMC"
noramp_params.optim_list=["k_0", "omega"]
num_params=5

param_ranges={
"E_0_range":[0.2, 0.3, 0.4, 0.5, 0.6],#np.linspace(estart+(de/2), ereverse-(de/2), num_params),
"k_0_range":[0.2,2, 20, 200, 2000],
"Ru_range":[1, 10, 100, 1000, 5000],#[10,100,500,1000,1500,2000,2500,3000],
"Cdl_range":[0,1e-7, 1e-6, 1e-5, 1e-4],
"CdlE1_range":[-2, -0.5,0, 0.5, 2],
"CdlE2_range":[-0.05, -0.01,0, 0.01, 0.05],
"CdlE3_range":[-0.01, -0.005,0, 0.005, 0.01],
"alpha_range":[0.1, 0.2,0.5, 0.7, 0.9],
}
param_keys=param_ranges.keys()
true_params=['E_0', 'k_0', 'Ru','Cdl', 'CdlE1','CdlE2','CdlE3','alpha']
start=time.time()
end=np.where(noramp_params.time_vec>0.01)
noramp_params.time_vec=np.array(noramp_params.time_vec)
cdl_count=0
voltages=noramp_params.define_voltages()
voltages=voltages
sort_idx=np.argsort(voltages)
voltages=np.multiply(voltages, noramp_params.nd_param.c_E0)
voltages=voltages[end]
num_times=5
#voltages=voltages[sort_idx]
bounds_vals=[200]
#plt.plot(voltage_results, current_results)
#plt.show()
Ru_range=param_ranges["Ru_range"]
debugging=False
surfaces=False
init_range=[-0.1, 0.1]
interval=abs(init_range[1]-init_range[0])/8
orig=init_range[0]
if surfaces ==True:
    time_idx=np.linspace(0, len(noramp_params.time_vec), num_times, endpoint=False)
else:
    time_idx=[1]
init_vals=np.linspace(-10, 10, 8)
print init_vals
means=[0.217740555023939, 357.3609131669447, 382.2055036116924, 5.440677575193328e-05, 0, 1.9319598326712751e-10, 8.94098189688349, 1.406746302896052, 0.9000000005558953]

list1=['E_0', 'k_0', 'Ru','Cdl', 'CdlE1','gamma','omega', 'phase','alpha']
noramp_params.optim_list=list1
time_series=noramp_params.test_vals(means, "timeseries", test=False)
plt.subplot(1,2,1)
plt.plot(voltage_results, current_results)
plt.subplot(1,2,2)
plt.plot(voltages, time_series[end])
plt.show()
params_1=np.zeros(len(list1))
for i in range(0, 8): #len(param_keys)
    for j in range(0, len(true_params)):
        if true_params[j] in param_keys[i]:
            if (true_params[j]=="Cdl") and (cdl_count>0):
                continue
            elif (true_params[j]=="Cdl"):
                cdl_count+=1
            parameter_name=true_params[j].strip()
            break

    if debugging ==True:
        for k in range(0, len(time_idx)):#num_params

            noramp_params.numerical_method="inverted"
            if surfaces ==True:
                    plt.subplot(2,4, i+1)
                    noramp_params.debug_time=noramp_params.time_vec[int(time_idx[k])]
                    noramp_params.initial_val=init_vals[i]
                    time_series=noramp_params.test_vals([3000, 8.94],"timeseries", test=False)
            #noramp_params.numerical_method=solver


            #
            elif surfaces==False:
                    for q in range(0, len(bounds_vals)):
                        plt.subplot(2,4, i+1)
                        noramp_params.bounds_val=bounds_vals[q]
                        time_series=noramp_params.simulate([Ru_range[i], 8.94], noramp_params.time_vec, "no", "timeseries", "no")
                        time_series=np.multiply(time_series, noramp_params.nd_param.c_I0)
                        plt.plot(noramp_params.time_vec, time_series, label=str(bounds_vals[q]))
                        plt.xlabel("Nondim current")
                        plt.ylabel("F(current)")
            else:
                init_times=np.linspace(orig, orig+interval, 5)
                for q in range(0, len(init_times)):
                    plt.subplot(2,4, i+1)
                    noramp_params.initial_val=init_times[q]
                    time_series=noramp_params.simulate([3000, 8.94], noramp_params.time_vec, "no", "timeseries", "no")
                    plt.plot(noramp_params.time_vec, time_series)
                orig=orig+interval

            plt.title("Ru="+str(Ru_range[i]))
    else:
        noramp_params.optim_list=[parameter_name, "omega"]
        position=true_params.index(parameter_name)
        plt.subplot(2,4,position+1)
        plt.title(parameter_name)
        #plt.plot(voltage_results, current_results)
        print noramp_params.simulation_options["numerical_method"]
        for k in range(0, num_params):
            for q in range(0, len(list1)):
                params_1[q]=noramp_params.dim_dict[list1[q]]
            print list(params_1)
            parameter_val=param_ranges[param_keys[i]][k]
            time_series=noramp_params.test_vals([parameter_val, 8.94], "timeseries", test=False)
            time_series=time_series[end]
            time_series=np.multiply(time_series, noramp_params.nd_param.c_I0)
            plt.plot(voltages, time_series, label=parameter_val)

        plt.legend()
        #noramp_params.nd_param.non_dimensionalise(parameter_name, param_list[parameter_name])

    noramp_params=single_electron(param_list, simulation_options, other_values)
plt.show()
