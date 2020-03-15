import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import sys
import matplotlib.ticker as ticker
import math
import scipy.stats as stat
dir_path = os.path.dirname(os.path.realpath(__file__))
slash_idx=[i for i in range(len(dir_path)) if dir_path[i]=="/"]
one_above=dir_path[:slash_idx[-1]]
sys.path.insert(1, one_above)
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
def det_subplots(value):
    if np.floor(np.sqrt(value))**2==value:
        return int(np.sqrt(value)), int(np.sqrt(value))
    if value<=10:
        start_val=2
    else:
        start_val=3

    rows=range(start_val, int(np.ceil(value/start_val)))
    for i in range(0, 10):
        modulos=np.array([value%x for x in rows])
        idx_0=(np.where(modulos==0))
        if len(idx_0[0])!=0:
            return int(rows[idx_0[0][-1]]), int(value/rows[idx_0[0][-1]])
        value+=1


def find(name, path, Electrode, skiprows=0):
    for root, dirs, files in os.walk(path):
        if Electrode in dirs:
            files=os.listdir(root+"/"+Electrode)
            if name in files:
                return np.loadtxt(root+"/"+Electrode+"/"+name, skiprows=skiprows)
def change_param(params, optim_list, parameter, value):
    param_list=copy.deepcopy(params)
    param_list[optim_list.index(parameter)]=value
    return param_list
def chain_appender(chains, param):
    new_chain=chains[0, :, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, :, param])
    return new_chain
def plot_params(titles, set_chain, positions=None, label=None, row=None, col=None):
    if positions==None:
        positions=range(0, len(titles))
    if row is None:
        row, col=det_subplots(len(titles))
    for i in range(0, len(titles)):
        axes=plt.subplot(row,col,i+1)
        plot_chain=chain_appender(set_chain, positions[i])
        if abs(np.mean(plot_chain))<0.001:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        elif abs(np.mean(plot_chain))<10:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        else:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        axes.hist(plot_chain, alpha=0.4,bins=20, stacked=True, label=label)
        axes.legend()
        lb, ub = axes.get_xlim( )
        axes.set_xticks(np.linspace(lb, ub, 4))
        axes.set_xlabel(titles[i])
        axes.set_ylabel('frequency')
        axes.set_title(graph_titles[i])
    plt.subplots_adjust(left=0.06, bottom=0.05, right=0.96, top=0.95, wspace=0.31, hspace=0.32)
def trace_plots(titles, chains, names, rhat=False, burn_in_thresh=0):
    row, col=det_subplots(len(titles))
    if rhat==True:
        rhat_vals=pints.rhat_all_params(chains[:, burn_in_thresh:, :])
    for i in range(0, len(titles)):
        axes=plt.subplot(row,col, i+1)
        for j in range(0, len(chains)):
            axes.plot(chains[j, :, i], label="Chain "+str(j), alpha=0.7)
        if abs(np.mean(chains[j, :, i]))<0.001:
            axes.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
        else:
            axes.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        if i==0:
            axes.legend()
        #lb, ub = axes.get_xlim( )
        #axes.set_xticks(np.linspace(lb, ub, 5))
        if rhat == True:
            axes.set_title(names[i]+ " (Rhat="+str(round(rhat_vals[i],3))+")")
        else:
            axes.set_title(names[i])
        axes.set_ylabel(titles[i])
        axes.set_xlabel('Iteration')
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.95, wspace=0.62, hspace=0.35)
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
    "cap_phase":"",
    'phase' : "",
    "alpha_mean": "",
    "alpha_std": "",
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
    "cap_phase":"$C_{dl}$ phase",
    "alpha_mean": "$\\alpha\\mu$",
    "alpha_std": "$\\alpha\\sigma$",
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
    "E_0":"Midpoint potential",
    "cap_phase":"Capacitance phase",
    "alpha_mean": "Symmetry factor mean",
    "alpha_std": "Symmetry factor standard deviation",
    'phase' : "Phase",
    "":"Experiment",
    "noise":"Noise",
}
#f=open(filename, "r")
params=["E0_mean", "E0_std", "k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha_std", "noise"]
params1=["E0_mean", "E0_std", "k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha_mean", "alpha_std", "noise"]
params=["E0_mean", "E0_std", "k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase","alpha_std", "noise"]
optim_list=params
titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
graph_titles=[Titles[x] for x in optim_list]
folder="MCMC_runs/omega_nondim/"
electrode="Yellow"
path=("/").join([dir_path , electrode, folder])
files=os.listdir(path)#
#fig=plt.figure(num=None, figsize=(12,9), dpi=120, facecolor='w', edgecolor='k')
concs1=["1e-1", "5e-1", "1e0", "15e-1", "3e0"]
concs=["__{0}__".format(x) for x in concs1]
concs=[str(x) for x in range(1, 6)]

vals=[0.1, 0.5, 1, 1.5, 3]
extension="_cv_high_ru_MCMC.run3"
desired_file="Noramp_"
file_numbers=["8_94","114", "209"]
#concs=["__{0}__".format(x) for x in file_numbers]
#positions_23=[params1.index(x) for x in params1]
#positions_24=[params1.index(x) for x in params2]
ns=[str(x) for x in range(1, 11)]
means=np.zeros((len(ns), len(params)))
stds=np.zeros((len(ns), len(params)))
nums=["_{0}_cv".format(x) for x in ns]
#for i in range(0, len(concs)):
runs=["run24"]
row, col=det_subplots(len(params1))
idx=-1
plt.rcParams.update({'font.size': 16})
def standardise(x, mean, standard, sqrt_N):
    return (x-mean)/(standard/sqrt_N)
normed_chains=np.zeros((len(titles), 50000*3, len(ns)))
for num in nums:
    idx+=1
    for filename in files:#
        run_check=[string in filename for string in runs]
        if True in run_check and num in filename:# and  True  in [x in filename for x in include]:
            #print(filename)
            number=filename[7]
            label_idx=run_check.index(True)
            chains=np.load(("/").join([electrode, folder, filename]))
            """
            if number=="2":
                plot_params(titles, chains[:, 25000:45000, :], positions=positions, label="Scan "+num)
            else:
                plot_params(titles, chains[:, 30000:, :], positions=positions, label="Scan "+num)
            """

            means[idx, :]=[np.mean(chain_appender(chains[:, 50000:, :], x)) for x in range(0, len(titles))]
            stds[idx, :]=[np.std(chain_appender(chains[:, 50000:, :], x)) for x in range(0, len(titles))]


global_mean=[np.mean(means[:, x]) for x in range(0, len(params))]
global_stds=[np.std(means[:, x]) for x in range(0, len(params))]

xlabel=ns

farad_params=["E0_mean", "E0_std", "k_0", "alpha_std"]
Cdl_params=["Cdl","CdlE1", "CdlE2"]
input_params=["omega","cap_phase","phase"]
misc_params=["Ru", "gamma", "noise"]
fig, axes=plt.subplots(2,2)
for i in range(0, len(titles)):
    subbed=np.subtract(means[: i],global_mean[i])
    subbed=[means[x, i]-global_mean[i] for x in range(0, len(ns))]
    normed_means=np.divide(subbed, global_stds[i]/np.sqrt(len(titles)))
    mean_ratio=np.divide(means[:,i], normed_means)
    scaled_error=np.divide(stds[:,i], mean_ratio)


    if params[i] in farad_params:
        ax=axes[0,0]
        ax.set_ylabel("Z-score")
    elif params[i] in Cdl_params:
        ax=axes[0,1]
    elif params[i] in input_params:
        ax=axes[1,0]
        ax.set_ylabel("Z-score")
        ax.set_xlabel("Experiment")
    elif params[i] in misc_params:
        ax=axes[1,1]
        ax.set_xlabel("Experiment")
    ax.plot(ns, normed_means, label=titles[i])
    ax.legend()
    fig.set_size_inches((6.5, 9))

plt.show()
save_path="Z-scores"fig.savefig(save_path, dpi=200)


            #print([np.mean(chains[1, 50000:, x]) for x in range(0, len(titles))])


            #trace_plots(titles, chains[:, :, :], graph_titles, rhat=True )
            #plt.show()
    #plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92, wspace=0.30, hspace=0.33)
#plot_params(titles, chains2)
