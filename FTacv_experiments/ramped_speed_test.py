import numpy as np
import matplotlib.pyplot as plt
import pints
from single_e_class_unified  import single_electron
from harmonics_plotter import harmonics
import pints.plot
import os
from multiplotter import multiplot
import math
dir_path = os.path.dirname(os.path.realpath(__file__))
results_dict="Inferred_params"
Electrode="Yellow"
run="Run_2"
counter=1
ramped_file="Ramped_3_cv_high_ru.ts"
ramp_data_class=single_electron(("/").join([dir_path, results_dict,Electrode,run, ramped_file]), {}, {}, {}, {}, False)
ramped_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase","phase", "alpha"]
ramp_data_class.def_optim_list(ramped_optim_list)
ramp_data_class.dim_dict["v_nondim"]=True
ramp_data_class.dim_dict["CdlE3"]=0
ramp_data_class.dim_dict["phase"]=0
ramp_data_class.dim_dict["cap_phase"]=0

#ramp_data_class.dim_dict["omega"]=8.95
ramp_data_class.dim_dict["original_omega"]=8.94
ramp_data_class.simulation_options["dispersion_bins"]=16
ramped_data_path=("/").join([dir_path, "experiment_data_2", "Experimental-120919", "Ramped", "Yellow"])
start_harm=2
end_harm=7
ramp_data_harm_class=harmonics(list(range(start_harm, end_harm)), ramp_data_class.nd_param.omega, 0.05)
harmonic_range=list(range(start_harm, end_harm))
ramp_data_class.harmonic_range=harmonic_range
def likelihood_func(harmonic_range, times, time_series, omega):
    filter_val=0.1
    frequency=np.fft.fftfreq(len(times), times[1]-times[0])
    hanning=np.hanning(len(times))
    time_series=np.multiply(time_series, hanning)
    desired_idx=np.where((frequency>((harmonic_range[0]*omega)-(omega*filter_val)))&(frequency<(harmonic_range[-1]*omega)+(omega*filter_val)))
    Y=(np.fft.fft(time_series))
    return Y[desired_idx],frequency[desired_idx]
method="timeseries"
counter=0
j=0
counter2=0
funcs=[np.real, np.imag]
labels=["Real", "Imaginary"]

i="2"
ramped_current_results=np.loadtxt(ramped_data_path+"/Yellow_Electrode_Ramped_"+str(i)+"_cv_current")
ramped_voltage_results=np.loadtxt(ramped_data_path+"/Yellow_Electrode_Ramped_"+str(i)+"_cv_voltage")[:,1]

ramped_time_results=ramped_current_results[:,0]
ramped_current_results=ramped_current_results[:,1]
ramp_data_class.other_values["experiment_time"]=ramped_time_results
ramp_data_class.other_values["experiment_current"]=ramped_current_results
ramp_data_class.other_values["experiment_voltage"]=ramped_voltage_results
ramp_fit=single_electron(None, dim_paramater_dictionary=ramp_data_class.dim_dict, simulation_options=ramp_data_class.simulation_options, other_values=ramp_data_class.other_values, param_bounds=ramp_data_class.param_bounds)
ramped_data_harmonics=ramp_data_harm_class.generate_harmonics(ramped_time_results, ramped_current_results)
file="Noramp_"+str(i)+"_cv_high_ru.run3_2"
master_optim_list=["E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","phase", "alpha"]
ramp_fit.def_optim_list(master_optim_list)
noramp_results=single_electron(("/").join([dir_path, results_dict,Electrode,"Run_3", file]))
noramp_harm=harmonics(list(range(start_harm, end_harm)), noramp_results.nd_param.omega, 0.05)
param_vals=([noramp_results.save_dict["params"][0][noramp_results.save_dict["optim_list"].index(key)] if  (key in noramp_results.save_dict["optim_list"]) else noramp_results.dim_dict[key] for key in master_optim_list])
noramp_results.def_optim_list(master_optim_list)

#print(param_vals)
cmaes_time_series=noramp_results.i_nondim(noramp_results.test_vals(param_vals, method))
current_results=noramp_results.i_nondim(noramp_results.other_values["experiment_current"])
voltage_results=noramp_results.e_nondim(noramp_results.other_values["experiment_voltage"])
time_results=noramp_results.t_nondim(noramp_results.other_values["experiment_time"])
ramp_fit.dim_dict["sampling_freq"]=1/200.0
cmaes_times=time_results#noramp_results.t_nondim(noramp_results.time_vec[noramp_results.time_idx:])
test_voltages=voltage_results
filtered_exp, filtered_exp_freq=likelihood_func(harmonic_range, cmaes_times, cmaes_time_series, noramp_results.nd_param.omega)
filtered_data, filtered_data_freq=likelihood_func(harmonic_range, time_results, current_results,  noramp_results.nd_param.omega)

#test_voltages=np.interp(cmaes_times, time_results, voltage_results)
ramped_param_vals=[noramp_results.dim_dict[x] for x in ramped_optim_list]
ramped_param_vals=[0.230288480903605996, 0.041226544633855209, 122.16602260400302, 742.0285824092344, 8.212777750912338e-05, 0.00295011584929738, -0.0005139373308230138, 7.330227126430568e-11, 8.885659818458482,0,0.5146469208293426]
values=[[0.23794865088573297, 0.012233237509352422, 127.31889407388738, 166.5594838874783, 7.735179022904297e-05, 0.002100658721047255, -0.0003313772993110888, 7.92266323442856e-11, 8.884799587792013, 0, 0.5999996893552537],
        [0.23382687987164236, 0.021002377079024245, 65.09191639664698, 670.8662113597322, 7.639003487570702e-05, 0.005475992161506471, -0.0005413285385403921, 7.597878417461771e-11, 8.884799587792013, 0, 0.55],
        [0.23422376456952138, 0.03260990314447788, 127.84406477439944, 866.5594838874783, 7.6536077935438e-05, 0.0028147360128534457, -0.00040069155469661145, 7.510752378483546e-11, 8.884799587792013,  0, 0.5999999355676289] ,
        [0.23388480903605996, 0.030956721349239304, 125.32474893610775, 866.5594838874783, 7.633158635035894e-05, 0.0025282530103625106, -0.00038943293495081674, 7.525190180668605e-11, 8.884799587792013, 0, 0.5952709809282939] ,
        [0.23251409481227692, 0.03738294965025301, 131.6521158525072, 866.5594838874783, 7.621231601962323e-05, 0.0026921147914714116, -0.00041118732365409173, 7.377445613940177e-11, 8.884799587792013,  0, 0.5999999998930867]]
values=[[0.22318477825106026, 0.05017555057827148, 299.99999708633436, 435.3085348834923, 7.656969863908582e-05, 0.0049162085192067315, -0.0004940984281294643, 7.132894746570913e-11, 8.884799587792013, 0, 0.46] ,
        [0.22370166148102177, 0.048594513570175274, 299.999999937947, 461.06579388366265, 7.668183755836e-05, 0.004462686615568902, -0.0004731551691502609, 7.195182931026769e-11, 8.884799587792013, 0, 0.48000000000000004] ,
        [0.23176101084840808, 0.030347974151306224, 107.53613561318956, 790.443764676077, 7.6820398998218e-05, 0.004474624678728965, -0.0004807322727249999, 7.466976561611421e-11, 8.884799587792013,0, 0.5] ,
        [0.23900539779698587, 0.012495118289935866, 56.188310062404945, 828.5655234790614, 7.657793838284919e-05, 0.005997891678714437, -0.0005672650186174608, 7.57075288471139e-11, 8.884799587792013,0, 0.52]]
values=[
        [0.23895989026998632, 0.010302865192620036, 60.4060892206524, 853.5788290922803, 7.688160708685767e-05, 0.004892910406443837, -0.0005177765609662311, 7.656559966657177e-11, 8.884799587792013, 0, 0.54] ,
        [0.236488634758585, 0.015094276283200609, 80.29560428116115,881.3101581551538, 7.734794864089517e-05, 0.0030106294613832726, -0.0004260879423058535, 7.746352664803931e-11, 8.884799587792013, 0, 0.56] ,
        [0.23390430905699255, 0.020859204824098317, 104.15659670702253, 884.6010926039512, 7.765425228548337e-05, 0.0016876381150878542, -0.0003603144768100069, 7.769052465288134e-11, 8.884799587792013,0, 0.5800000000000001] ,
        [0.2298173759746179, 0.0303946808562666, 141.1395680666635, 932.148247550758, 7.773297133386808e-05, 0.0011426547005235627, -0.00033131533359411804, 7.637796505553124e-11, 8.884799587792013, 0, 0.6100000000000001] ,
        [0.23110722113853555, 0.027345924242786514, 128.71431443672168, 854.1702606093729, 7.77430224842735e-05, 0.0011777514485454693, -0.0003337620915544407, 7.694184659356877e-11, 8.884799587792013, 0, 0.6000000000000001] ,
]
#values=[0.24326919974213232, 106.44472632014336, 856.2055228333725, 7.730195292672856e-05, 0.0021848907053882236, -0.00038393545481927183, 8.118068324955082e-11, 8.884799587792013,0,0.577740128733653]



#axes=multiplot(1, 1, **{"harmonic_position":0, "num_harmonics":5, "orientation":"landscape", "plot_width":5})
#[0.23363378517047495, 0.03481010462713212, 125.2418680556816, 866.5594838874783, 7.694171200331618e-05, 0.003209611999861764, -0.0004263185805571494, 7.476933579891946e-11, 8.884799587792013, 5.02832295276801, 0.5999996422725197] ,
#ramped_param_vals=[0.2146905833892083, 0.03054987719700828, 97.73281867896537, 607.3604383735801, 9.855333151006176e-05, 0.002360092744665423, -0.0006167247956683124, 2.816480836302093e-10, 8.884799587792013, 3.7667007049769547, 5.036133769048488, 0.6175763030220238]
#ramped_param_vals=[0.2105491969962405, 0.060492147323936776, 109.9494203694055, 667.8257252271911, 9.034055522879727e-05, 0.002655104265433323, -0.0005653310639016844, 1.1167514700671097e-10, 8.884680131969208, 1.232169931719665, 0.46318222889196947]
#ramped_param_vals=[0.228308296004341, 0.04678928704199179, 109.9494215030595, 753.3526020630202, 7.39149998133382e-05, 0.0032451274206835596, -0.00046254359814023406, 8.063249835588253e-11, 8.884745411351478, 2.236534676329661, 0.4631822303721044]
#ramped_param_vals=[0.22285754799270494, 0.04694179204711807, 247.23685999878182, 512.1240571770394, 7.706083559507504e-05, 0.0028478498239308467, -0.00041016659783854256, 7.257026464642658e-11, 8.885659818458482,0, 0.6999999975394792]
ramp_fit.harmonic_range=harmonic_range
#axes_keys=sorted(axes.axes_dict.keys())
for i in range(0, 10):
    ramped_time_series=ramp_fit.i_nondim(ramp_fit.test_vals(values[0], "timeseries"))

    """
    alpha_val=round(ramp_fit.dim_dict["alpha"],3)
    ramped_times=ramp_fit.t_nondim(ramp_fit.time_vec[ramp_fit.time_idx:])
    ramped_harmonics=ramp_data_harm_class.generate_harmonics(ramped_times, ramped_time_series)
    results=np.loadtxt(dir_path+"/experiment_data_2/Experimental-120919/Ramped/Yellow/"+"Yellow_Electrode_Ramped_2_cv_current")
    time_results=results[:,0]
    current_results=results[:,1]
    sinusoid_data_harmonics=noramp_harm.generate_harmonics(time_results, current_results)
    sinusoid_harmonics=noramp_harm.generate_harmonics(cmaes_times, cmaes_time_series)
    ramp_fit.simulation_options["method"]="dcv"
    dcv_plot=ramp_fit.e_nondim(ramp_fit.define_voltages())
    ramp_fit.simulation_options["method"]="ramped"


    for harm_counter in range(0, len(ramped_harmonics)):
        print(j)
        ax=axes.axes_dict["row1"][j]
        #print(j, j%((len(values)//2)*len(ramped_harmonics)))
        ax.plot(ramped_times, (ramped_harmonics[harm_counter,:]*1e6), label="Simultation")
        ax.plot(ramped_time_results, (ramped_data_harmonics[harm_counter,:]*1e6), label="Data", alpha=0.5)
        ax2=ax.twinx()
        ax2.set_ylabel(harmonic_range[harm_counter], rotation=0)
        ax2.set_yticks([])
        if harm_counter==0:
            ax.set_title("No dispersion")
            #plt.axvline(ramped_param_vals[0], color="black", linestyle="--", label="$E_0$")
        if harm_counter%3==2:
            ax.set_ylabel("Current($\mu$A)")
        if harm_counter==(len(ramped_harmonics)-1):
            ax.set_xlabel("Time(s)")
            ax.legend()
        else:
            ax.set_xticks([])
        #sinusoidal_axes[j].plot(cmaes_times, (sinusoid_harmonics[harm_counter,:]*1e6), label="Simulation")
        #sinusoidal_axes[j].plot(time_results, (sinusoid_data_harmonics[harm_counter,:]*1e6), label="Data")
        j+=1

"""
print(len(ramped_time_series), len(ramped_time_results))
print(np.mean(ramp_fit.time_array), np.std(ramp_fit.time_array))
plt.show()
