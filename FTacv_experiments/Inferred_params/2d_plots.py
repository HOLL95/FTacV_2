import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import scipy.stats as stats
from matplotlib.ticker import FormatStrFormatter
from PIL import Image
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
def chain_appender(chains, param):
    print(len(chains))
    if len(chains)>20:
        return chains[:, param]
    new_chain=chains[0, :, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, :, param])
    return new_chain
all_params=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "cap_phase","phase", "alpha_std", "noise"]
#all_params=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "phase","cap_phase", "alpha"]
optim_list=['E0_mean', "E0_std",'k_0',"Ru", "phase", "cap_phase","alpha_std"]
all_params=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "cap_phase","phase", "alpha_mean","alpha_std", "noise"]
#all_params=['E0_mean', "E0_std",'k_0',"Ru","Cdl", "CdlE1", "CdlE2",'gamma',"omega", "phase","cap_phase", "alpha"]
optim_list=optim_list
positions=[all_params.index(x) for x in optim_list]
#positions[-1]=positions[-1]-1

titles=[fancy_names[x]+"("+unit_dict[x]+")" if (unit_dict[x]!="") else fancy_names[x] for x in optim_list]
n_param=len(titles)
file="Noramp_2_cv_high_ru.run3_2_MCMC_run15"
file="Noramp_2_cv_high_ru_alpha_0.61_MCMC_run15"
file="Noramp_2_cv_high_ru_MCMC.run24"
file="Noramp_2_cv_high_ru.run3_2_MCMC_run21"
file="Noramp_2_cv_high_ru_MCMC.run23"
file="Noramp_2_cv_high_ru.run3_2_MCMC_run7"
file="Noramp_2_cv_high_ru_MCMC.run23"
folder="MCMC_runs/omega_nondim/"
electrode="Yellow"
path=("/").join([dir_path , electrode, folder])
files=os.listdir(path)

#chain_result[1, :, ]=chain_result[2, :, :]
#chain_result=chain_result[:2, :, :]
#chains=chains[1, 10000:, :]
#chains=chains[:, 40000:, :]

def plot_kde_1d(x, ax, num=None):
    """ Creates a 1d histogram and an estimate of the PDF using KDE. """
    xmin = np.min(x)
    xmax = np.max(x)
    ax.set_xlim(xmin, xmax)
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
plt.rcParams.update({'font.size': 16})
def plot_kde_2d(x, y, ax):
    # Get minimum and maximum values
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # Plot values
    #values = np.vstack([x, y])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    ax.scatter(x, y, s=0.5, alpha=0.5)#, cmap=plt.cm.Blues, extent=[xmin, xmax, ymin, ymax])

    # Create grid
    #xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    #positions = np.vstack([xx.ravel(), yy.ravel()])

    # Get kernel density estimate and plot contours
    #kernel = stats.gaussian_kde(values)
    #f = np.reshape(kernel(positions).T, xx.shape)
    #ax.contourf(xx, yy, f, cmap='Blues')
    #ax.contour(xx, yy, f, colors='k')
    ax.locator_params(nbins=2)
    #if np.std(x)<4e-4:
    #    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    #    ax.set_xticks([xmin, np.mean(x)])
    #if np.std(y)<4e-4:
    #    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    # Fix aspect ratio
    #print(((xmax - xmin)/ (ymax - ymin)))
    #ax.set_aspect(0.25*(xmax - xmin)/ (ymax - ymin))
run="run24"
num=["_"+str(x)+"_" for x in range(10, 11)]
for number in num:
    for file in files:
        if number in file and run in file:
            print(number)
            print(file)
            fig_size=(12,12)
            fig, ax=plt.subplots(n_param, n_param)
            chain_result=np.load(("/").join([dir_path, "Yellow", "MCMC_runs","omega_nondim",file ]))
            for q in range(0, len(chain_result)):#len(chain_result)
                #chains=chain_result[:, 80000:, :]
                chain_len=len(chain_result[q, :, 0])
                chains=chain_result[q, 50000:, :]
                for i in range(0,n_param):
                    for j in range(0, n_param):
                        if i==j:
                            axes=ax[i,j]


                            #axes.set_yticks([])
                            ax1=axes.twinx()
                            #plot_kde_1d(chain_appender(chains, positions[j]), ax=axes)
                            plot_kde_1d(chains[:,positions[j]], ax=axes, num=("Chain "+str(q+1)))
                            if i==0:
                                axes.legend(loc="center left", bbox_to_anchor=(1.75, 0.5))
                            ticks=axes.get_yticks()
                            #labels=axes.get_yticklabes([])
                            axes.set_yticks([])
                            ax1.set_yticks(ticks)
                            if q==0:
                                #ax1.set_yticks(ticks)
                                ax1.set_ylabel("frequency")

                            #x1.set_ytick_labels()
                            #if i==0:
                            #    axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                            #axes.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))

                            #axes.set_yticks([])
                        elif i<j:
                            ax[i,j].axis('off')
                        else:
                            axes=ax[i,j]
                            #plot_kde_2d(chain_appender(chains, positions[j]), chain_appender(chains, positions[i]), ax=axes)
                            plot_kde_2d(chains[:,positions[j]],chains[:,positions[i]], ax=axes)

                        if i!=0:
                            ax[i, 0].set_ylabel(titles[i])
                        if i<n_param-1:
                            ax[i,j].set_xticklabels([])#
                        if j>0 and i!=j:
                            ax[i,j].set_yticklabels([])
                        if j!=n_param:
                            ax[-1, i].set_xlabel(titles[i])
                            plt.setp( ax[-1, i].xaxis.get_majorticklabels(), rotation=15 )

            plt.subplots_adjust(left=0.1, bottom=0.1, right=0.91, top=0.98, wspace=0.2, hspace=0.12)
            fig = plt.gcf()
            fig.set_size_inches((14,9))
            plt.show()
            save_path="Alice_"+str(number)+"_2d_plots.png"
            fig.savefig(save_path, dpi=500)
            img = Image.open(save_path)
            basewidth = float(img.size[0])//2
            wpercent = (basewidth/float(img.size[0]))
            hsize = int((float(img.size[1])*float(wpercent)))
            img = img.resize((int(basewidth),hsize), Image.ANTIALIAS)
            img.save(save_path, "PNG", quality=95, dpi=(500, 500))
