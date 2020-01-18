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
run="Run_5"
concs=["1e-1M", "1e0M"]
file_numbers=[str(x) for x in range(1, 4)]
#figure=multiplot(4, 4, **{"harmonic_position":3, "num_harmonics":5, "orientation":"portrait", "fourier_position":2, "plot_width":5})
#keys=sorted(figure.axes_dict.keys())
plot_counter=0
h_counter=0
f_counter=0
other_files=["6", "9"]
def RMSE(series1, series2):
    return np.sqrt((np.sum(1/(len(series1))*np.power(np.subtract(series1, series2),2))))
for i in np.arange(0.45, 0.68, 0.01):

    file="Noramp_2_cv_high_ru_alpha_"+str(round(i, 3))

    method="timeseries"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))

    #For 1e0M no. 1
    #[0.4437541435808381, 0.04672244339040179, 169.82717561197694, 1499.9999999884658, 3.326734620433954e-05, -0.0024550352309892637, 0.003444105740388075, 1.6931489825942052e-10, 8.94186753976662, 4.218129657177546, 5.5734366297761]
    #For 1e0M no.2
    #[0.4480443313256469, 0.04578316505621488, 186.54661732132573, 1499.999998740905, 3.2749457208458154e-05, -0.0037426290649702418, 0.0033541225901685782, 1.5457625945960093e-10, 8.941954712163895, 4.2264148667766035, 5.54918407688431]
    #For 1e0M no. 3
    #[0.44902174254951466, 0.04649482119604907, 167.77163066745334, 1499.9999899319125, 3.244942640679048e-05, -0.0015943375724840197, 0.003256602506039867, 1.5805338795893364e-10, 8.941715397492652, 4.227414643240183, 5.562298818419109]


    noramp_results.simulation_options["dispersion_bins"]=32
    noramp_results.param_bounds["alpha"]=[0.1, 0.9]
    param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
    noramp_results.def_optim_list(master_optim_list)
    cmaes_time=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))
    current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])#[0::dec_amount]
    voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])#[0::dec_amount]
    print_vals=np.append(param_vals, RMSE(current_results, cmaes_time)*1e6)
    print(list(print_vals), ",")
    #plt.plot(voltage_results, current_results)
    #plt.plot(voltage_results, cmaes_time, alpha=0.7)
    #plt.show()

    #param_vals=[0.448042594478853, 0.04578302993986619, 186.55619822678182, 1499.9999952103126, 3.274964218833196e-05, -0.0037437099369155846, 0.0033541585110203383, 1.5457010853658904e-10, 8.941956053700354, 4.226410309725663, 5.549175826705732]

    #noramp_results.dim_dict["alpha"]=0.5
    #noramp_results.simulation_options["label"]="MCMC"
    #param_vals[2]=150

    #param_vals=[0.17989588529462708, 0.0561105224748124, 180.3580615548242, 1266.0203230196744, 1.8968981153440156e-05, 0.0194709122755512, 0.0001216629532319758, 1.5349480973755286e-10, 209.77882631191702, 4.551912859108585, 6.12365070126766, 0.4499006211167625]


plt.show()
