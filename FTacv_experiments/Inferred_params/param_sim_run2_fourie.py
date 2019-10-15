import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle
dir_path = os.path.dirname(os.path.realpath(__file__))
slash_idx=[i for i in range(len(dir_path)) if dir_path[i]=="/"]
one_above=dir_path[:slash_idx[-1]]
sys.path.insert(1, one_above)
from single_e_class_unified import single_electron
from harmonics_plotter import harmonics
Electrode="Yellow"
path=("/").join([dir_path , Electrode])
files=os.listdir(path)#
file_numbers=[]
counter=1
cols=["low"]
rows=[str(x) for x in range(1,4)]
#fig, ax=plt.subplots(len(rows),len(cols))
optim_list=['E0_mean', "E0_std","k_0","Ru", "Cdl","CdlE1","CdlE2", "gamma", "omega","phase", "cap_phase", "alpha"]
def RMSE(vec1, vec2):
    return np.sqrt((np.sum(np.power(np.subtract(vec1, vec2),2)))/len(vec1))
results_array=[]
name_list=[]
counter=0
for lcv_1 in range(0, len(cols)):
    for lcv_2 in range(0, len(rows)):
        #dot_idx=names.index(".")
        filename=("_").join(["Noramp", rows[lcv_2], "cv", cols[lcv_1], "ru"])+".fourier1"

        #file_numbers.append(names[names.index(".")-1])
        result=single_electron(path+"/"+filename)
        result.simulation_options["dispersion_bins"]=20
        voltage_results=result.other_values["experiment_voltage"]
        current_results=result.other_values["experiment_current"]
        time_results=result.other_values["experiment_time"]
        harm_class=harmonics(range(3, 9,1), result.dim_dict["omega"], 0.1)
        likelihood_func=result.kaiser_filter(current_results)
        #for i in range(0, len(self.save_dict[params])):
        #if lcv_2==0:
        #    ax[lcv_2, lcv_1].set_title(cols[lcv_1]+ " Resistance ")
        #if lcv_2==len(rows)-1:
        #    ax[lcv_2, lcv_1].set_xlabel("Nondim voltage")
        #if lcv_1==0:
        #    ax[lcv_2, lcv_1].set_ylabel(rows[lcv_2])
        result.simulation_options["likelihood"]="timeseries"
        ru_idx=result.save_dict["optim_list"].index("Ru")
        #result.save_dict["params"][0][ru_idx]=65
        simulations=(result.saved_param_simulate(result.save_dict["params"][0]))
        print_list=[result.dim_dict[key] for key in optim_list]
        #print_list=np.append(print_list, RMSE(simulations,result.other_values["experiment_current"])*result.nd_param.c_I0)
        results_array.append(list(print_list))
        name_list.append(cols[lcv_1].capitalize()+"\\ "+ "$R_u$"+"\\ " +rows[lcv_2])
        counter+=1#
        data_harms=harm_class.generate_harmonics(time_results, current_results)
        sim_harms=harm_class.generate_harmonics(time_results, simulations)
        harm_class.plot_harmonics(time_results, "phased", data=data_harms, sim=sim_harms)
        #plt.plot(simulations)
        #plt.plot(likelihood_func, alpha=0.7)
        #plt.show()
        #inv_func=harm_class.inv_objective_fun(result.kaiser_filter, simulations)
        #inv_func_data=harm_class.inv_objective_fun(result.kaiser_filter, likelihood_func)
        #plt.subplot(len(cols), len(rows), counter)
        #plt.plot(np.fft.ifft(simulations))
        #plt.plot(np.fft.ifft(likelihood_func))
        #plt.plot(np.subtract(np.fft.ifft(simulations),np.fft.ifft(likelihood_func)))
        #plt.plot(result.other_values["experiment_voltage"], simulations)
        #plt.plot(result.other_values["experiment_voltage"], result.other_values["experiment_current"], alpha=0.7)
        #plt.plot(result.other_values["experiment_voltage"], np.subtract(simulations,result.other_values["experiment_current"]))
#fig.text(0.05, 0.5, 'Nondim current', ha='center', va='center', rotation='vertical')
plt.show()
for i in range(0, len(results_array)):
    print results_array[i]
#with open("alice_yellow_params.pkl", "wb") as f:
#    pickle.dump([name_list, results_array], f)
