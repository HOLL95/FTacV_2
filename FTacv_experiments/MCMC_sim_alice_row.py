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
run="Run_6"
concs=["1e-1M", "1e0M"]
file_numbers=[str(x) for x in range(1, 4)]
num_h=7
figure=multiplot(1, 3, **{"harmonic_position":2, "num_harmonics":num_h, "orientation":"portrait", "fourier_position":1, "plot_width":5})
#plt.show()
keys=sorted(figure.axes_dict.keys())
plot_counter=0
h_counter=0
f_counter=0
other_files=["6", "9"]
def RMSE(series1, series2):
    return np.sqrt((np.sum(1/(len(series1))*np.power(np.subtract(series1, series2),2))))
for i in range(2,3):
    if str(i) in other_files:
        file="Noramp_"+str(i)+"_cv_high_ru.run3_4"
    else:
        file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
    file="Noramp_"+str(i)+"_cv_high_ru_alpha_disp"
    method="timeseries"
    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
    noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,run, file]))

    #For 1e0M no. 1
    #[0.4437541435808381, 0.04672244339040179, 169.82717561197694, 1499.9999999884658, 3.326734620433954e-05, -0.0024550352309892637, 0.003444105740388075, 1.6931489825942052e-10, 8.94186753976662, 4.218129657177546, 5.5734366297761]
    #For 1e0M no.2
    #[0.4480443313256469, 0.04578316505621488, 186.54661732132573, 1499.999998740905, 3.2749457208458154e-05, -0.0037426290649702418, 0.0033541225901685782, 1.5457625945960093e-10, 8.941954712163895, 4.2264148667766035, 5.54918407688431]
    #For 1e0M no. 3
    #[0.44902174254951466, 0.04649482119604907, 167.77163066745334, 1499.9999899319125, 3.244942640679048e-05, -0.0015943375724840197, 0.003256602506039867, 1.5805338795893364e-10, 8.941715397492652, 4.227414643240183, 5.562298818419109]


    noramp_results.simulation_options["dispersion_bins"]=16
    noramp_results.simulation_options["GH_quadrature"]=True
    noramp_results.simulation_options["alpha_dispersion"]=True
    noramp_results.param_bounds["alpha"]=[0.58, 0.65]
    param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])

    master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha_mean", "alpha_std"]
    param_vals=[0.2437136033540455, 0.041916959869871935, 136.91071718538336, 535.9775743351688, 7.719208829248059e-05, 0.002419672072987142, -0.00039798693040548304, 7.389022703687543e-11, 8.940524733361801, 4.357762581087894, 4.974566671527844, 0.6065682978303066, 0.17229158014845414]
    #param_vals=[0.24704003369743605, 127.22888060724983, 899.0593759401338, 7.759349521072281e-05, 0.0018792925278760739, -0.0003331168610431201, 8.105946935824678e-11, 8.940643709990447, 4.457735184257261, 5.135033708944046, 0.596434662870841]
    noramp_results.def_optim_list(master_optim_list)
    cmaes_time=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))
    current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])#[0::dec_amount]
    voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])#[0::dec_amount]
    #plt.plot(voltage_results, current_results)
    #plt.plot(voltage_results, cmaes_time, alpha=0.7)
    #plt.show()
    print_vals=np.append(param_vals, RMSE(current_results, cmaes_time)*1e6)
    print(list(print_vals), ",")
    #param_vals=[0.448042594478853, 0.04578302993986619, 186.55619822678182, 1499.9999952103126, 3.274964218833196e-05, -0.0037437099369155846, 0.0033541585110203383, 1.5457010853658904e-10, 8.941956053700354, 4.226410309725663, 5.549175826705732]

    #noramp_results.dim_dict["alpha"]=0.5
    #noramp_results.simulation_options["label"]="MCMC"
    #param_vals[2]=150

    #param_vals=[0.17989588529462708, 0.0561105224748124, 180.3580615548242, 1266.0203230196744, 1.8968981153440156e-05, 0.0194709122755512, 0.0001216629532319758, 1.5349480973755286e-10, 209.77882631191702, 4.551912859108585, 6.12365070126766, 0.4499006211167625]


    harms=harmonics(range(1, 1+num_h),noramp_results.dim_dict["omega"] , 0.05)
    #noramp_results.simulation_options["likelihood"]=method
    #noramp_results.simulation_options["dispersion_bins"]=16
    plt.rcParams.update({'axes.labelsize': 16})
    mpl.rcParams['axes.labelsize'] = 12
    dec_amount=16
    print(len(figure.axes_dict[keys[0]]))
    time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])#[0::dec_amount]
    figure.axes_dict[keys[0]][0].plot(voltage_results*1e3, (cmaes_time)*1e3, label="Sim")
    figure.axes_dict[keys[0]][0].plot(voltage_results*1e3,(current_results)*1e3, alpha=0.5, label="Data")
    figure.axes_dict[keys[0]][0].set_xlabel("Voltage(mV)")
    figure.axes_dict[keys[0]][0].set_ylabel("Current(mA)")
    figure.axes_dict[keys[0]][0].plot(voltage_results*1e3, np.subtract(cmaes_time*1e3, current_results*1e3), label="Residual")
    figure.axes_dict[keys[0]][0].legend()
    plot_counter+=1
    harms.harmonic_selecter(figure.axes_dict[keys[1]][0], cmaes_time, time_results,  box=False, arg=np.real, line_label="Sim")
    harms.harmonic_selecter(figure.axes_dict[keys[1]][0],current_results, time_results,  box=False, arg=np.real, line_label="Data")
    figure.axes_dict[keys[1]][0].set_xlabel("Frequency(Hz)")
    figure.axes_dict[keys[1]][0].set_ylabel("Real")
    figure.axes_dict[keys[1]][0].legend(loc="lower left")
    f_counter+=1
    harms.harmonic_selecter(figure.axes_dict[keys[1]][1], cmaes_time, time_results, box=False, arg=np.imag,line_label="Sim")
    harms.harmonic_selecter(figure.axes_dict[keys[1]][1],  current_results,  time_results,box=False, arg=np.imag, line_label="Data")
    figure.axes_dict[keys[1]][1].set_xlabel("Frequency(Hz)")
    figure.axes_dict[keys[1]][1].set_ylabel("Imaginary")
    #figure.axes_dict[keys[0]][f_counter].legend(loc=2)
    f_counter+=1
    data_harms=harms.generate_harmonics(time_results, cmaes_time)
    exp_harms=harms.generate_harmonics(time_results, current_results)
    for q in range(0, len(data_harms)):
        figure.axes_dict[keys[2]][h_counter].plot(voltage_results*1e3, data_harms[q]*1e6, label="Sim")
        figure.axes_dict[keys[2]][h_counter].plot(voltage_results*1e3, exp_harms[q]*1e6, alpha=0.5, label="Data")
        lb, ub = figure.axes_dict[keys[2]][h_counter].get_ylim( )
        figure.axes_dict[keys[2]][h_counter].set_yticks([round(x,1) for x in np.linspace(0.5*lb, 0.5*ub, 2)])
        ax2=figure.axes_dict[keys[2]][h_counter].twinx()
        ax2.set_ylabel(harms.harmonics[q], rotation=0)
        ax2.set_yticks([])
        if q%len(data_harms)==2:
            figure.axes_dict[keys[2]][h_counter].set_ylabel("Current($\mu$A)")
        if q==len(data_harms)-1:
            figure.axes_dict[keys[2]][h_counter].set_xlabel("Voltage(mV)")
        #if q==0:
            figure.axes_dict[keys[2]][h_counter].legend(loc="lower right")
        h_counter+=1


plt.show()
