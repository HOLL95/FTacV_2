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
import pints.plot
from rdp import rdp
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
Experiment="Varying_freq"
folder="Noramp"

type="current"
type2="voltage"
rhat_mean=2
path=("/").join([dir_path, data_path, folder, Experiment])
files= os.listdir(path)
def binary_file_reader(filename):
    file=open(filename, "rb")
    binary_array=[]
    for line in file:
        string=line.decode("latin1")
        binary_array.append([float(string[0:string.index("\t")]), float(string[string.index("\t")+len("\t"):string.index("\r")])])
    return np.array(binary_array)
file_numbers=["8_94","114", "209"]
file_nos=["_{0}_".format(x) for x in file_numbers]
dec_amounts=[32, 1, 1]
omegas=[8.94]#, 114.14, 209.18]
for counter in range(0, len(file_nos)):
    desired_conc=file_nos[counter]
    binary_files={}
    for filename in files:
        for type in types:
            if type in filename and desired_conc in filename:
                print(filename)
                binary_files[type]=path+"/"+filename



    experiment_data={}
    for type in types:
        experiment_data[type]=binary_file_reader(binary_files[type])
    dec_amount=dec_amounts[counter]
    #rdp_results=rdp(list(zip(experiment_data["current"][0::dec_amount,0],experiment_data["current"][0::dec_amount, 1])))
    #rdp_results=np.array(rdp_results)
    #current_results=np.convolve(experiment_data["current"][:, 1], np.ones((dec_amount,))/dec_amount, mode="valid")
    #current_results1=current_results[0::dec_amount]
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
        'omega':omegas[counter],#8.88480830076,  #    (frequency Hz)
        "original_omega":omegas[counter],
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
        'phase' : 4.949075135495371,
        "time_end": None,
        'num_peaks': 11,
    }
    solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
    likelihood_options=["timeseries", "fourier"]
    time_start=1/(param_list["omega"])

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
        'Ru': [0, 1e4],  #     (uncompensated resistance ohms)
        'Cdl': [0,1e-4], #(capacitance parameters)
        'CdlE1': [-0.05,0.05],#0.000653657774506,
        'CdlE2': [-0.01,0.01],#0.000245772700637,
        'CdlE3': [-0.01,0.01],#1.10053945995e-06,
        'gamma': [1e-11,1e-9],
        'k_0': [10, 1e3], #(reaction rate s-1)
        'alpha': [0.4, 0.6],
        "cap_phase":[math.pi, 2*math.pi],
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
    blank_fit.def_optim_list(["Ru","Cdl","CdlE1", "CdlE2", "cap_phase", "omega"])
    blank_fit.nd_param.phase=(3*math.pi/2)+(math.pi/2)
    param_set_1=[2.7371615612353777e-07, 1.377386734451008e-05, -0.023521133369626402, 0.0007344175954793371, 4.281633502258979, 8.941765261346148]
    param_set_2=[999.9999983932938, 1.3523706647576712e-05, -0.021871525678528143, 0.0007168124240944443,3*math.pi/2, 8.941765261346148]
    sets=[param_set_1, param_set_2]
    labels=["Fitted phase", "Set phase"]
    plt.plot(blank_fit.e_nondim(voltage_results), blank_fit.i_nondim(current_results)*1e6, label="Data")
    for j in range(0, len(sets)):
        time_series=blank_fit.i_nondim(blank_fit.test_vals(sets[j], likelihood="timeseries"))*1e6
        plt.plot(blank_fit.e_nondim(voltage_results), time_series, linestyle="--", label=labels[j])
        plt.legend()
        plt.xlabel("Voltage(V)")
        plt.ylabel("Current($\\mu$A)")
    fig=plt.gcf()
    fig.set_size_inches((3.25, 4.5))
    plt.show()
    save_path="blank_phases.png"
    fig.savefig(save_path, dpi=500)

    """
    ax2=ax[0]
    ax2.plot(blank_fit.t_nondim(time_results), blank_fit.e_nondim(blank_fit.define_voltages(True)*1e3), color="red", label="Predicted capacitance phase+$\\frac{\pi}{2}$", linestyle="--")

    ax2.set_ylabel("Voltage(mV)")
    ax1=ax2.twinx()
    ax1.plot(blank_fit.t_nondim(time_results), blank_fit.i_nondim(current_results)*1e3)
    ax1.set_ylabel("Current(mA)")
    ax2.set_xlabel("Time(s)")
    legend_1 =ax2.legend(loc=1, borderaxespad=1.)
    legend_1.remove()
    legend_2=ax1.legend(loc=1, borderaxespad=1.)
    ax1.add_artist(legend_1)
    legend_2.remove()
    ax3=ax[1]
    ax3.set_ylabel("Voltage(mV)")

    blank_fit.nd_param.phase=4.3593306496513256+(math.pi/2)
    ax3.plot(blank_fit.t_nondim(time_results), blank_fit.e_nondim(blank_fit.define_voltages(True)*1e3), color="red", label="Fitted capacitance phase+$\\frac{\pi}{2}$", linestyle="--")

    ax4=ax[1].twinx()
    ax4.plot(blank_fit.t_nondim(time_results), blank_fit.i_nondim(current_results)*1e3)
    ax4.set_ylabel("Current(mA)")

    legend_1 =ax3.legend(loc=1, borderaxespad=1.)
    legend_1.remove()
    legend_2=ax4.legend(loc=1, borderaxespad=1.)
    ax4.add_artist(legend_1)
    legend_2.remove()
    ax3.set_xlabel("Time(s)")
    plt.show()
    plt.plot(time_results, voltage_results)
    #plt.plot(time_results, blank_fit.define_voltages()[blank_fit.time_idx:])
    #ax=plt.gca()
    #ax2=ax.twinx()
    #ax2.plot(time_results, current_results, color="red")
    plt.show()
    """
    blank_fit.dim_dict["cap_phase"]=3*math.pi/2
    blank_fit.dim_dict["omega"]=8.941765261346148
    blank_fit.def_optim_list(["Ru","Cdl","CdlE1", "CdlE2", "cap_phase"])
    cmaes_problem=pints.SingleOutputProblem(blank_fit, time_results, current_results)
    score = pints.SumOfSquaresError(cmaes_problem)#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
    print(list([np.zeros(len(blank_fit.optim_list))]), list([np.ones(len(blank_fit.optim_list))]))
    CMAES_boundaries=pints.RectangularBoundaries(list([np.zeros(len(blank_fit.optim_list))]), list([np.ones(len(blank_fit.optim_list))]))
    blank_fit.simulation_options["label"]="cmaes"
    x0=abs(np.random.rand(blank_fit.n_parameters()))#blank_fit.change_norm_group(gc4_3_low_ru, "norm")

    cmaes_fitting=pints.OptimisationController(score, x0, sigma0=None, boundaries=CMAES_boundaries, method=pints.CMAES)
    cmaes_fitting.set_max_unchanged_iterations(iterations=200, threshold=1e-3)
    cmaes_fitting.set_parallel(True)
    found_parameters, found_value=cmaes_fitting.run()
    cmaes_results=blank_fit.change_norm_group(found_parameters, "un_norm")
    print(list(cmaes_results))
    cmaes_time=blank_fit.test_vals(cmaes_results, likelihood="timeseries")
    plt.plot(voltage_results, cmaes_time)
    plt.plot(voltage_results, current_results)
    plt.show()
