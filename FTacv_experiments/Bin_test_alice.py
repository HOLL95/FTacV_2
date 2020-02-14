import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pints
from single_e_class_unified  import single_electron
from scipy.stats import norm
import pints.plot
import os
import math
from multiplotter import multiplot
from harmonics_plotter import harmonics
import time
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_6"
concs=["1e-1M", "1e0M"]
file_numbers=[str(x) for x in range(1, 4)]
#figure=multiplot(4, 4, **{"harmonic_position":3, "num_harmonics":5, "orientation":"portrait", "fourier_position":2, "plot_width":5})
#keys=sorted(figure.axes_dict.keys())
plot_counter=0
h_counter=0
f_counter=0
other_files=["6", "9"]
def normal_vals(mean, std):
    
def RMSE(series1, series2):
    return np.sqrt((np.sum(1/(len(series1))*np.power(np.subtract(series1, series2),2))))
for i in range(2, 11):

    file="Noramp_"+str(i)+"_cv_high_ru_alpha_disp"

    method="timeseries"


    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))
    bin_range=[5, 16, 50]
    error=np.zeros(len(bin_range))
    for j in range(0, len(bin_range)):
        noramp_results.simulation_options["dispersion_bins"]=bin_range[j]
        #For 1e0M no. 1
        #[0.4437541435808381, 0.04672244339040179, 169.82717561197694, 1499.9999999884658, 3.326734620433954e-05, -0.0024550352309892637, 0.003444105740388075, 1.6931489825942052e-10, 8.94186753976662, 4.218129657177546, 5.5734366297761]
        #For 1e0M no.2
        #[0.4480443313256469, 0.04578316505621488, 186.54661732132573, 1499.999998740905, 3.2749457208458154e-05, -0.0037426290649702418, 0.0033541225901685782, 1.5457625945960093e-10, 8.941954712163895, 4.2264148667766035, 5.54918407688431]
        #For 1e0M no. 3
        #[0.44902174254951466, 0.04649482119604907, 167.77163066745334, 1499.9999899319125, 3.244942640679048e-05, -0.0015943375724840197, 0.003256602506039867, 1.5805338795893364e-10, 8.941715397492652, 4.227414643240183, 5.562298818419109]

        master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha_mean", "alpha_std"]
        #noramp_results.param_bounds["alpha"]=[0.1, 0.9]
        param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
        #master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha_mean", "alpha_std"]
        #master_optim_list=["E_0","k0_shape", "k0_scale","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase", "phase", "alpha"]
        #noramp_results.dim_dict["alpha_mean"]=0.5
        #noramp_results.dim_dict["alpha_std"]=0.1
        #noramp_results.param_bounds["alpha_mean"]=[0.4, 0.6]
        #noramp_results.param_bounds["alpha_std"]=[1e-3, 0.3]

        #noramp_results.simulation_options["GH_quadrature"]=True

        noramp_results.def_optim_list(master_optim_list)
        #param_vals=[0.23424750529310384, 0.020320742816477166, 97.0885510867113, 667.9290934936241, 7.776160856516875e-05, 0.0010637389431854757, -0.0003320205908225078, 7.769167120665223e-11, 8.940519592519038, 4.418204867343626, 5.0732166228450275, 0.5500000001067766, 0.10672951494126706]
        #param_vals=[0.242821641458727, 0.5796871925249325, 49.95306048243109, 709.4025928368603, 7.702176462584796e-05, 0.004000430062620887, -0.0004796319064888713, 7.706597103131586e-11, 8.940522305966446, 4.405392753944093, 5.134525768028457, 0.5576193566666604]
        #param_vals=[0.23546569801184455, 0.04061615022925082, 140.42709391306067, 539.1331187810706, 7.71808303441534e-05, 0.0025956881531416925, -0.0004060837908083821, 7.375023058466999e-11, 8.940524602855044, 4.3586160431434555, 4.974285257652122, 0.5715519870619699, 0.15509617967563671]
        #param_vals=[0.24312173517024177, 0.036978313503345725, 126.83248590634965, 598.883974612922, 7.702632337023351e-05, 0.0034403582922662435, -0.0004433702323469006, 7.45077819233449e-11, 8.940524543980159, 4.37500416257147, 5.013410834767164, 0.5999999996254317]

        #param_vals=values=[0.236212535397608, 0.07718482329772401, 3.3878479936851233, 4.241119820972833e-07, 2.8617144117465527e-05, 0.1499999999841155, -0.007969322490141095, 3.865452240417434e-11, 8.940544410322712, 4.38668688522284, 3.9538725896604654, 0.5277937191659627]
        start=time.time()
        cmaes_time=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))
        #print(time.time()-start)
        current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])#[0::dec_amount]
        voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])#[0::dec_amount]
        #print_vals=np.append(param_vals, RMSE(current_results, cmaes_time)*1e6)

        #print(list(print_vals), ",")
        #
        plt.plot(voltage_results, cmaes_time, alpha=0.7, label=str(bin_range[j])+" Bins")
    plt.xlabel("Voltage(V)")
    plt.ylabel()
    #plt.plot(voltage_results, current_results, alpha=0.5, label="Data")
    plt.legend()

    plt.show()

    #param_vals=[0.448042594478853, 0.04578302993986619, 186.55619822678182, 1499.9999952103126, 3.274964218833196e-05, -0.0037437099369155846, 0.0033541585110203383, 1.5457010853658904e-10, 8.941956053700354, 4.226410309725663, 5.549175826705732]

    #noramp_results.dim_dict["alpha"]=0.5
    #noramp_results.simulation_options["label"]="MCMC"
    #param_vals[2]=150

    #param_vals=[0.17989588529462708, 0.0561105224748124, 180.3580615548242, 1266.0203230196744, 1.8968981153440156e-05, 0.0194709122755512, 0.0001216629532319758, 1.5349480973755286e-10, 209.77882631191702, 4.551912859108585, 6.12365070126766, 0.4499006211167625]


plt.show()
