import numpy as np
import matplotlib.pyplot as plt
import math
import time
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
import os
import pickle
import pints
dir_path = os.path.dirname(os.path.realpath(__file__))
types=["current", "voltage"]
exp="Experimental-120919"
bla="Blank-131119"

exp_type=bla
if exp_type==bla:
    extra="Blank/"
else:
    extra=""
data_path="experiment_data_2/"+exp_type
Experiment="Varying_salt"
folder="Noramp"


type="current"
type2="voltage"

path=("/").join([dir_path, data_path, folder, Experiment])
files= os.listdir(path)
desired_conc="1"
binary_files={}
for filename in files:
    for type in types:
        if type in filename and desired_conc in filename:
            binary_files[type]=path+"/"+filename


def binary_file_reader(filename):
    file=open(filename, "rb")
    binary_array=[]
    for line in file:
        string=line.decode("latin1")
        binary_array.append([float(string[0:string.index("\t")]), float(string[string.index("\t")+len("\t"):string.index("\r")])])
    return np.array(binary_array)
experiment_data={}
for type in types:
    experiment_data[type]=binary_file_reader(binary_files[type])
dec_amount=32
current_results1=experiment_data["current"][0::dec_amount, 1]
time_results1=experiment_data["current"][0::dec_amount,0]
voltage_results1=experiment_data["voltage"][0::dec_amount,1]
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    "E_0":0.2,
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-6, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 0,
    "original_gamma":1e-10,        # (surface coverage per unit area)
    'k_0': 10, #(reaction rate s-1)
    'alpha': 0.5,
    "E0_mean":0.2,
    "E0_std": 0.09,
    "k0_shape":0.954,
    "k0_loc":100,
    "k0_scale":50,
    "k0_range":1e3,
    "cap_phase":0,
    'sampling_freq' : (1.0/400),
    'phase' : 3*(math.pi/2),
    "time_end": None,
    'num_peaks': 50
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
time_start=2/(param_list["omega"])

simulation_options={
    "no_transient":time_start,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":16,
    "test": False,
    "method": "sinusoidal",
    "phase_only":False,
    "likelihood":likelihood_options[0],
    "numerical_method": solver_list[1],
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":list(range(3,9,1)),
    "experiment_time": time_results1,
    "experiment_current": current_results1,
    "experiment_voltage":voltage_results1,
    "bounds_val":20,
}
param_bounds={
    'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.95*param_list['omega'],1.05*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 700],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-0.05,0.15],#0.000653657774506,
    'CdlE2': [-0.01,0.01],#0.000245772700637,
    'CdlE3': [-0.01,0.01],#1.10053945995e-06,
    'gamma': [1e-11,1e-9],
    'k_0': [10, 1e3], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    "cap_phase":[0, 2*math.pi],
    "E0_mean":[0.19, 0.208],
    "E0_std": [0.001, 0.2],
    "k0_shape":[0,2],
    "k0_loc":[0, 1e3],
    "k0_scale":[0,1e3],
    "k0_range":[1e2, 1e4],
    'phase' : [0, 2*math.pi]
}
blank_fit=single_electron(None, param_list, simulation_options, other_values, param_bounds)
blank_fit.define_boundaries(param_bounds)
time_results=blank_fit.other_values["experiment_time"]
current_results=blank_fit.other_values["experiment_current"]
voltage_results=blank_fit.other_values["experiment_voltage"]
voltages=blank_fit.define_voltages()
blank_fit.def_optim_list(["Ru","Cdl","CdlE1","CdlE2","cap_phase"])
inferred_params=[2.7381994349086846e-08, 2.05933286880881e-05, -0.013538849600882703, 0.00023437171907022408, 4.3593306496513256]
plt.plot(voltage_results, current_results)
for i in [0, 10, 200, 300]:
    inferred_params[0]=i
    time_series=blank_fit.test_vals(inferred_params, "timeseries")
    plt.plot(voltage_results, time_series, alpha=0.5, label=i)
plt.legend()
plt.show()

cmaes_problem=pints.SingleOutputProblem(blank_fit, time_results, current_results)
score = pints.SumOfSquaresError(cmaes_problem)#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
print(list([np.zeros(len(blank_fit.optim_list))]), list([np.ones(len(blank_fit.optim_list))]))
CMAES_boundaries=pints.RectangularBoundaries(list([np.zeros(len(blank_fit.optim_list))]), list([np.ones(len(blank_fit.optim_list))]))
blank_fit.simulation_options["label"]="cmaes"
x0=abs(np.random.rand(blank_fit.n_parameters()))#blank_fit.change_norm_group(gc4_3_low_ru, "norm")
cmaes_fitting=pints.OptimisationController(score, x0, sigma0=None, boundaries=CMAES_boundaries, method=pints.CMAES)
cmaes_fitting.set_max_unchanged_iterations(iterations=200, threshold=1e-3)
cmaes_fitting.set_parallel(False)
found_parameters, found_value=cmaes_fitting.run()
cmaes_results=blank_fit.change_norm_group(found_parameters, "un_norm")
print(list(cmaes_results))
cmaes_time=blank_fit.test_vals(cmaes_results, likelihood="timeseries", test=True)
plt.plot(voltage_results, current_results)
plt.plot(voltage_results, cmaes_time)
plt.show()
