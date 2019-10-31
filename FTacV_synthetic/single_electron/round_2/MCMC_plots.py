import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import sys
import matplotlib.ticker as ticker
import math
dir_path = os.path.dirname(os.path.realpath(__file__))
slash_idx=[i for i in range(len(dir_path)) if dir_path[i]=="/"]
one_above=dir_path[:slash_idx[-1]]
sys.path.insert(1, one_above)



def change_param(params, optim_list, parameter, value):
    param_list=copy.deepcopy(params)
    param_list[optim_list.index(parameter)]=value
    return param_list
def chain_appender(chains, param):
    new_chain=chains[0, :, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, :, param])
    return new_chain
def plot_params(titles, set_chain, param_indexes, labels, boundaries):
    for i in range(0, len(titles)):
        print(i)
        axes=plt.subplot(2,4,i+1)
        plot_chain=chain_appender(set_chain, param_indexes[i])#set_chain[:, param_indexes[i]]#
        if abs(np.mean(plot_chain))<0.001:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3e'))
        else:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        axes.hist(plot_chain, alpha=0.4,bins=20, stacked=True, edgecolor='black', label=labels)
        #axes.set_xlim(boundaries[0][i], boundaries[1][i])
        lb, ub = axes.get_xlim( )

        axes.set_xticks(np.linspace(lb, ub, 4))
        axes.set_xlabel(titles[i])
        axes.set_ylabel('frequency')
        axes.set_title(graph_titles[i])
        axes.legend()


unit_dict={
    "E_0": "V",
    'E_start': "V", #(starting dc voltage - V)
    'E_reverse': "V",
    'omega':"Hz",#8.88480830076,  #    (frequency Hz)
    'd_E': "V",   #(ac voltage amplitude - V) freq_range[j],#
    'v': '$s^{-1}$',   #       (scan rate s^-1)
    'area': '$cm^{2}$', #(electrode surface area cm^2)
    'Ru': "$\\Omega$",  #     (uncompensated resistance ohms)
    'Cdl': "F", #(capacitance parameters)
    'CdlE1': "",#0.000653657774506,
    'CdlE2': "",#0.000245772700637,
    'CdlE3': "",#1.10053945995e-06,
    'gamma': 'mol cm^{-2}$',
    'k_0': 's^{-1}$', #(reaction rate s-1)
    'alpha': "",
    "E0_mean":"V",
    "E0_std": "V",
    "k0_shape":"",
    "k0_loc":"",
    "k0_scale":"",
    "cap_phase":"rads",
    'phase' : "rads",
    "":"",
    "noise":"",
}
fancy_names={
    "E_0": '$E^0$',
    'E_start': '$E_{start}$', #(starting dc voltage - V)
    'E_reverse': '$E_{reverse}$',
    'omega':'$\\omega$',#8.88480830076,  #    (frequency Hz)
    'd_E': "$\\Delta E$",   #(ac voltage amplitude - V) freq_range[j],#
    'v': "v",   #       (scan rate s^-1)
    'area': "Area", #(electrode surface area cm^2)
    'Ru': "Ru",  #     (uncompensated resistance ohms)
    'Cdl': "$C_{dl}$", #(capacitance parameters)
    'CdlE1': "$C_{dlE1}$",#0.000653657774506,
    'CdlE2': "$C_{dlE2}$",#0.000245772700637,
    'CdlE3': "$C_{dlE3}$",#1.10053945995e-06,
    'gamma': '$\\Gamma',
    'k_0': '$k_0', #(reaction rate s-1)
    'alpha': "$\\alpha$",
    "E0_mean":"$E^0 \\mu$",
    "E0_std": "$E^0 \\sigma$",
    "cap_phase":"Capacitance phase",
    'phase' : "Phase",
    "":"Experiment",
    "noise":"$\sigma$",
}
Titles={
    'omega':'Input frequency',#8.88480830076,  #    (frequency Hz)
    'd_E': "Amplitude",   #(ac voltage amplitude - V) freq_range[j],#
    'v': "Scan rate",   #       (scan rate s^-1)
    'area': "Area", #(electrode surface area cm^2)
    'Ru': "Uncompensated resistance",  #     (uncompensated resistance ohms)
    'Cdl': "Linear capacitance", #(capacitance parameters)
    'CdlE1': "First order capacitance",#0.000653657774506,
    'CdlE2': "Second order capacitance",#0.000245772700637,
    'CdlE3': "Third order capacitance",#1.10053945995e-06,
    'gamma': 'Surface coverage',
    'k_0': 'Rate constant', #(reaction rate s-1)
    'alpha': "Symmetry factor",
    "E0_mean":"Themodynamic mean",
    "E0_std": "Thermodynamic standard deviation",
    "E_0": "Midpoint potential",
    "cap_phase":"Capacitance phase",
    'phase' : "Phase",
    "":"Experiment",
    "noise":"Noise",
}
#f=open(filename, "r")
all_params=["E_0", "k_0", "Ru", "gamma", "Cdl", "CdlE1", "CdlE2", "phase", "cap_phase", "alpha", "noise"]
optim_list=["E_0", "k_0", "Ru", "gamma", "Cdl", "phase", "cap_phase", "alpha"]
param_idx=[all_params.index(x) for x in optim_list]
print(param_idx)
true_params=[0.25, 200, 300, 1e-10, 1e-5, -1e-4, 1e-4, 3*math.pi/2, 3*math.pi/2, 0.5]

desired_vals=[true_params[x] for x in param_idx]
boundaries=[np.multiply(desired_vals, 0.99), np.multiply(desired_vals, 1.01)]
boundaries=np.sort(boundaries, axis=0)
titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
graph_titles=[Titles[x] for x in optim_list]



j=0
#fig=plt.figure(num=None, figsize=(12,9), dpi=120, facecolor='w', edgecolor='k')
filename_list=[]
number_list=[]
num_params=len(optim_list)
electrode="Yellow"
omegas=[10, 100]
omega_label=["$\omega$="+str(x) for x in omegas]
for j in range(0, len(param_idx)):
    plt.subplot(2,4, j+1)
    plt.axvline(true_params[param_idx[j]], linestyle="--", color="black")


for i in range(0, len(omegas)):
    filename="Synthetic_MCMC_"+str(omegas[i])+"_run2"
    chains=np.load(filename)
    #pints.plot.trace(chains)
    plot_params(titles, chains[:, 5000:, :], param_idx, omega_label[i], boundaries)
    print(pints.rhat_all_params(chains[:, 5000:, :]))
plt.show()
