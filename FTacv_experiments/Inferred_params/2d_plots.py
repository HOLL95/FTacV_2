import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import scipy.stats as stats
from matplotlib.ticker import FormatStrFormatter
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
def chain_appender(chains, param):
    new_chain=chains[0, 30000:, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, 30000:, param])
    return new_chain
all_params=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "phase","cap_phase"]
optim_list=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "phase","cap_phase"]
positions=[all_params.index(x) for x in optim_list]
#positions[-1]=positions[-1]-1
plt.rcParams.update({'font.size': 12})
titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
n_param=len(titles)
chains=np.load(("/").join([dir_path, "Yellow", "MCMC_runs","Noramp_3_timeseriesMCMC_1_run1"]))
chains=chains[:, :, :]
fig_size=(12,12)
fig, ax=plt.subplots(n_param, n_param)
def plot_kde_1d(x, ax):
    """ Creates a 1d histogram and an estimate of the PDF using KDE. """
    xmin = np.min(x)
    xmax = np.max(x)
    x1 = np.linspace(xmin, xmax, 100)
    x2 = np.linspace(xmin, xmax, 50)
    ax.hist(x, bins=x2)
    if np.mean(x)<0.001:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    elif np.mean(x)<1:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.4f'))
    else:
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_xticks(np.linspace(xmin, xmax, 4))
def plot_kde_2d(x, y, ax):
    # Get minimum and maximum values
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # Plot values
    values = np.vstack([x, y])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(np.linspace(xmin, xmax, 3))
    ax.imshow(np.rot90(values), cmap=plt.cm.Blues, extent=[xmin, xmax, ymin, ymax])

    # Create grid
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])

    # Get kernel density estimate and plot contours
    kernel = stats.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    ax.contourf(xx, yy, f, cmap='Blues')
    ax.contour(xx, yy, f, colors='k')
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
    # Fix aspect ratio
    print(((xmax - xmin)/ (ymax - ymin)))
    ax.set_aspect((xmax - xmin)/ (ymax - ymin))
for i in range(0,n_param):
    for j in range(0, n_param):
        if i==j:
            axes=ax[i,j]
            if titles[i]=='Cdl':
                axes.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
            plot_kde_1d(chain_appender(chains, positions[j]), ax=axes)
            ax[i, j].set_ylabel("frequency")
        elif i<j:
            ax[i,j].axis('off')
        else:
            axes=ax[i,j]
            plot_kde_2d(chain_appender(chains, positions[j]), chain_appender(chains, positions[i]), ax=axes)
        if i!=0:
            ax[i, 0].set_ylabel(titles[i])
        if j!=n_param:
            ax[-1, i].set_xlabel(titles[i])

plt.show()
