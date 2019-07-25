import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import matplotlib.ticker as ticker
def chain_appender(chains, param):
    new_chain=chains[0, 5000:, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, 5000:, param])
    return new_chain
files=os.listdir('.')
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
    "E0_mean":"E^0 \\mu",
    "E0_std": "E^0 \\sigma",
    "cap_phase":"Capacitance phase",
    'phase' : "Phase",
    "":"Experiment",
    "noise":"$\sigma$",
}

#f=open(filename, "r")
optim_list=['E_0','k_0', 'Ru',"Cdl","CdlE1", "CdlE2",'gamma', "omega","phase", "cap_phase", "alpha", "noise"]
titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
true_params=[-0.4, 1e1, 0.000134*1e4,20.0]
noise_vals=[0.01, 0.02, 0.005, 0]
counter=-1
extension=".txt"
chain_select=1
j=0
#fig=plt.figure(num=None, figsize=(12,9), dpi=120, facecolor='w', edgecolor='k')
filename_list=[]
number_list=[]

chains=np.load('GC4_MCMC_chains_constrained.txt')
chains2=np.load("GC4_MCMC_chains_alpha_unconru.txt")
def plot_params(titles, chain):
    for i in range(0, len(titles)):
        print i
        axes=plt.subplot(3,4,i+1)
        plot_chain=chain_appender(chains, i)
        if abs(np.mean(plot_chain))<0.001:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3e'))
        else:
            axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        axes.hist(plot_chain, alpha=0.4,bins=20, stacked=True, edgecolor='black')
        lb, ub = axes.get_xlim( )
        axes.set_xticks(np.linspace(lb, ub, 4))
        axes.set_xlabel(titles[i])
        axes.set_ylabel('frequency')
    #if titles[i]!='$\sigma$':
    #    axes.axvline(true_params[i],color="black", linestyle="--")
    #axes.legend()
            #axes.set_xlim((0.9*true_params[i], 1.1*true_params[i]))
    #axes=plt.subplot(3,4,i+1)
    #axes.axis('off')
    plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92, wspace=0.30)
plot_params(titles,chains)
optim_list=['E_0','k_0', 'Ru',"Cdl","CdlE1", "CdlE2",'gamma', "omega","phase", "cap_phase", "alpha","noise"]
titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
plt.show()
for i in range(0, len(titles)):
    print titles[i]

    axes=plt.subplot(3,4,i+1)
    if optim_list[i]=="alpha":
        plt.axvline(0.5, color="black", linestyle="--")
        axes.set_xlabel(titles[i])
        axes.set_ylabel('frequency')
        continue
    plot_chain=chain_appender(chains2, j)
    if abs(np.mean(plot_chain))<0.001:
        axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3e'))
    else:
        axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    axes.hist(plot_chain, alpha=0.4,bins=20, stacked=True, edgecolor='black')
    lb, ub = axes.get_xlim( )
    axes.set_xticks(np.linspace(lb, ub, 4))
    axes.set_xlabel(titles[i])
    axes.set_ylabel('frequency')
    j+=1
plt.show()
