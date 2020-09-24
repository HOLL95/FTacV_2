import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import scipy.stats as stats
from matplotlib.ticker import FormatStrFormatter
import itertools
from mpl_toolkits import mplot3d
dir_path = os.path.dirname(os.path.realpath(__file__))
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
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
    "cap_phase":"Capacitance phase",
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
def chain_appender(chains, param):
    print(len(chains))
    if len(chains)>20:
        return chains[:, param]
    new_chain=chains[0, :, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, :, param])
    return new_chain
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
all_params=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "phase","cap_phase", "alpha_std"]
#all_params=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "phase","cap_phase", "alpha"]
optim_list=['E0_mean', "E0_std",'k_0',"Ru", "alpha_std"]
positions=[all_params.index(x) for x in optim_list]
#positions[-1]=positions[-1]-1
plt.rcParams.update({'font.size': 12})
titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
n_param=len(titles)
file="Noramp_2_cv_high_ru.run3_2_MCMC_run15"
file="Noramp_2_cv_high_ru_alpha_0.61_MCMC_run15"
file="Noramp_2_cv_high_ru_MCMC.run24"
file="Noramp_2_cv_high_ru.run3_2_MCMC_run21"
file="Noramp_2_cv_high_ru_MCMC.run23"
file="Noramp_2_cv_high_ru.run3_2_MCMC_run7"
file="Noramp_2_cv_high_ru_MCMC.run24"
chain_result=np.load(("/").join([dir_path, "Yellow", "MCMC_runs","omega_nondim",file ]))

#chain_result[1, :, ]=chain_result[2, :, :]
#chain_result=chain_result[:2, :, :]
#chains=chains[1, 10000:, :]
#chains=chains[:, 40000:, :]
#fig_size=(12,12)
#fig, ax=plt.subplots(n_param, n_param)
def plot_kde_1d(x, ax, num=None):
    """ Creates a 1d histogram and an estimate of the PDF using KDE. """
    xmin = np.min(x)
    xmax = np.max(x)
    x1 = np.linspace(xmin, xmax, 100)
    x2 = np.linspace(xmin, xmax, 50)
    ax.hist(x, bins=x2, label=num)
    """if np.mean(x)<0.001:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    elif np.mean(x)<1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    else:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))"""
    #ax.set_xticks(np.linspace(xmin, xmax, 4))
def plot_kde_2d(x, y, ax):
    # Get minimum and maximum values
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # Plot values
    #values = np.vstack([x, y])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    #ax.set_xticks(np.linspace(xmin, xmax, 3))
    ax.scatter(x, y, s=0.5, alpha=0.5)#, cmap=plt.cm.Blues, extent=[xmin, xmax, ymin, ymax])

    # Create grid
    #xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    #positions = np.vstack([xx.ravel(), yy.ravel()])

    # Get kernel density estimate and plot contours
    #kernel = stats.gaussian_kde(values)
    #f = np.reshape(kernel(positions).T, xx.shape)
    #ax.contourf(xx, yy, f, cmap='Blues')
    #ax.contour(xx, yy, f, colors='k')
    """
    if np.mean(x)<0.001:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    elif np.mean(x)<1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    else:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if np.mean(y)<0.001:
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    elif np.mean(y)<1:
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    else:
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    """
    # Fix aspect ratio
    #print(((xmax - xmin)/ (ymax - ymin)))
    #ax.set_aspect(0.25*(xmax - xmin)/ (ymax - ymin))
combinations=list(itertools.combinations(optim_list, r=3))
num_row, num_col=det_subplots(len(combinations))

#fig, ax=plt.subplots(num_row, num_col)
for m in range(0,len(combinations)):
    fig=plt.figure()
    ax=fig.add_subplot(1, 1, 1, projection="3d")
    positions=[all_params.index(x) for x in combinations[m]]
    chains=[chain_appender(chain_result[:, 50000:, :], x) for x in positions]
    ax.scatter(chains[0],chains[1], chains[2], s=0.5)
    ax.set_xlabel(fancy_names[combinations[m][0]]+unit_dict[combinations[m][0]])
    ax.set_ylabel(fancy_names[combinations[m][1]]+unit_dict[combinations[m][1]])
    ax.set_zlabel(fancy_names[combinations[m][2]]+unit_dict[combinations[m][2]])
    for angle in range(0, 360):
        ax.view_init(30, angle)
        plt.draw()
        plt.pause(.001)

#plt.subplots_adjust(left=0.07, bottom=0.08, right=0.96, top=0.98, wspace=0.44, hspace=0.32)