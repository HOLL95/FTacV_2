import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pints
from single_e_class_unified  import single_electron
import pints.plot
import os
import math
from multiplotter import multiplot
from harmonics_plotter import harmonics
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_3"
concs=["1e-1M", "1e0M"]
file_numbers=[str(x) for x in range(1, 4)]
plot_params=["E0_mean", "E0_std","k_0","Ru","gamma", "alpha"]
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase","alpha"]
num_harms=5
figure=multiplot(num_rows=2, num_cols=len(plot_params)/2, **{"harmonic_position":[0,1], "num_harmonics":num_harms, "orientation":"landscape",  "plot_width":5})


keys=sorted(figure.axes_dict.keys())

fancy_names={
    "E_0": '$E^0$',
    'E_start': '$E_{start}$', #(starting dc voltage - V)
    'E_reverse': '$E_{reverse}$',
    'omega':'$\\omega$',#8.88480830076,  #    (frequency Hz)
    'd_E': "$\\Delta E$",   #(ac voltage amplitude - V) freq_range[j],#
    'v': "v",   #       (scan rate s^-1)
    'area': "Area", #(electrode surface area cm^2)
    'Ru': "$R_{u}$",  #     (uncompensated resistance ohms)
    'Cdl': "$C_{dl}$", #(capacitance parameters)
    'CdlE1': "$C_{dlE1}$",#0.000653657774506,
    'CdlE2': "$C_{dlE2}$",#0.000245772700637,
    'CdlE3': "$C_{dlE3}$",#1.10053945995e-06,
    'gamma': '$\\Gamma$',
    'k_0': '$k_0$', #(reaction rate s-1)
    'alpha': "$\\alpha$",
    "E0_mean":"$E^0 \\mu$",
    "E0_std": "$E^0 \\sigma$",
    "cap_phase":"Capacitance phase",
    'phase' : "Phase",
    "":"Experiment",
    "noise":"$\sigma$",
}
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
    'gamma': '$mol\ cm^{-2}$',
    'k_0': '$s^{-1}$', #(reaction rate s-1)
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
plot_counter=0
h_counter=0
f_counter=0
other_files=["6", "9"]
def RMSE(series1, series2):
    return np.sqrt((np.sum(1/(len(series1))*np.power(np.subtract(series1, series2),2))))
param_bounds={
    'E_0':[0.3, 0.6],#[param_list['E_start'],param_list['E_reverse']],
    'Ru': [0,1000],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-0.05,0.15],#0.000653657774506,
    'CdlE2': [-0.01,0.01],#0.000245772700637,
    'CdlE3': [-0.01,0.01],#1.10053945995e-06,
    'gamma': [8e-11,1e-9],
    'k_0': [50, 200], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    "cap_phase":[5*math.pi/4, 2*math.pi],
    "E0_mean":[0.1, 0.4],
    "E0_std": [0.05, 0.2],
    'phase' : [5*math.pi/4, 2*math.pi]
}
i=2
if str(i) in other_files:
    file="Noramp_"+str(i)+"_cv_high_ru.run3_4"
else:
    file="Noramp_"+str(i)+"_cv_high_ru.run3_2"

method="timeseries"

noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))

#For 1e0M no. 1
#[0.4437541435808381, 0.04672244339040179, 169.82717561197694, 1499.9999999884658, 3.326734620433954e-05, -0.0024550352309892637, 0.003444105740388075, 1.6931489825942052e-10, 8.94186753976662, 4.218129657177546, 5.5734366297761]
#For 1e0M no.2
#[0.4480443313256469, 0.04578316505621488, 186.54661732132573, 1499.999998740905, 3.2749457208458154e-05, -0.0037426290649702418, 0.0033541225901685782, 1.5457625945960093e-10, 8.941954712163895, 4.2264148667766035, 5.54918407688431]
#For 1e0M no. 3
#[0.44902174254951466, 0.04649482119604907, 167.77163066745334, 1499.9999899319125, 3.244942640679048e-05, -0.0015943375724840197, 0.003256602506039867, 1.5805338795893364e-10, 8.941715397492652, 4.227414643240183, 5.562298818419109]


noramp_results.simulation_options["dispersion_bins"]=50
param_vals=[0.25, 0.1, 100, 100, 1e-10, 0.5]
noramp_results.dim_dict["Cdl"]=1e-5
noramp_results.dim_dict["CdlE1"]=0
noramp_results.dim_dict["CdlE2"]=0
noramp_results.def_optim_list(plot_params)
cmaes_time=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))
current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])
voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])
time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])

print_vals=np.append(param_vals, RMSE(current_results, cmaes_time)*1e6)
print(list(print_vals), ",")
#param_vals=[0.448042594478853, 0.04578302993986619, 186.55619822678182, 1499.9999952103126, 3.274964218833196e-05, -0.0037437099369155846, 0.0033541585110203383, 1.5457010853658904e-10, 8.941956053700354, 4.226410309725663, 5.549175826705732]

#noramp_results.dim_dict["alpha"]=0.5
#noramp_results.simulation_options["label"]="MCMC"
#param_vals[2]=150

#param_vals=[0.17989588529462708, 0.0561105224748124, 180.3580615548242, 1266.0203230196744, 1.8968981153440156e-05, 0.0194709122755512, 0.0001216629532319758, 1.5349480973755286e-10, 209.77882631191702, 4.551912859108585, 6.12365070126766, 0.4499006211167625]

start_harm=3
harms=harmonics(range(start_harm, start_harm+num_harms),noramp_results.dim_dict["omega"] , 0.05)
#noramp_results.simulation_options["likelihood"]=method
#noramp_results.simulation_options["dispersion_bins"]=16
num_scans=3
plt.rcParams.update({'axes.labelsize': 16})
mpl.rcParams['axes.labelsize'] = 12
for j in range(0, len(plot_params)):
    col_idx=int(j%(len(plot_params)/2))
    row_idx=int(np.floor(2*j/len(plot_params)))
    new_params=np.copy(param_vals)
    values=np.linspace(param_bounds[plot_params[j]][0], param_bounds[plot_params[j]][1],num_scans)

    for q in range(0, num_scans):
        new_params[j]=values[q]
        cmaes_time=noramp_results.i_nondim(noramp_results.test_vals(new_params, method))
        scan_harmonics=harms.generate_harmonics(time_results, cmaes_time)
        for z in range(0, num_harms):
            ax=figure.axes_dict[keys[row_idx]][(num_harms*col_idx)+z]
            if z==0:
                ax.set_title(fancy_names[plot_params[j]])
            if abs(values[q])>0.009:
                num_label=str(round(values[q], 3))
            else:
                num_label='%.2E' % values[q]
            ax2=ax.twinx()
            ax2.set_yticks([])
            ax2.set_ylabel(start_harm+z, rotation=0)
            ax.plot(voltage_results, scan_harmonics[z]*1e6, label=num_label+unit_dict[plot_params[j]])
            if z<num_harms-1:
                ax.set_xticks([])
            else:
                ax.set_xlabel("Voltage(V)")
            if z==(num_harms-1)/2:
                ax.set_ylabel("Current($\mu$A)")
                ax.legend(loc="center right")



    #for q in range(0, num_harms):




plt.show()
