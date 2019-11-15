import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_noramp  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
params_for_opt=[]

dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Black"
Method ="N_"
type="current"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
        print data
length_list=[2e4, 3e4, 4e4]
dec_list=[8,16, 32]
repeat_num=2
desired_length=int(length_list[0])
dec_amount=8
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
de=300e-3
estart=260e-3-de
ereverse=estart+2*de

param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 2.7,  #     (uncompensated resistance ohms)
    'Cdl': 0.00134, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-11,          # (surface coverage per unit area)
    'k_0': 10000.0, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
    'time_end':1000,
    'num_peaks': 50
}


de_novo=True
param_list['E_0']=(param_list['E_reverse']-param_list['E_start'])/2
harmonic_range=np.arange(4,9,1)
noramp_fit=single_electron(param_list, params_for_opt, harmonic_range, 0.5)
noramp_fit.filter_val=noramp_fit.nd_param.c_T0
noramp_fit.label="cmaes"
time_results=time_results[:desired_length]/noramp_fit.nd_param.c_T0#
current_results=current_results[:desired_length]/noramp_fit.nd_param.c_I0#
noramp_fit.throw_error=True
#plt.plot(time_results, current_results)
#plt.show()
noramp_fit.experimental_times=time_results
noramp_fit.time_vec=time_results#np.arange(0, noramp_fit.nd_param.sampling_freq*5000, noramp_fit.nd_param.sampling_freq)#
print noramp_fit.nd_param.sampling_freq
signal_length=len(current_results)
noramp_fit.num_points=signal_length
frequencies=np.fft.fftfreq(signal_length, noramp_fit.time_vec[1]-noramp_fit.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
noramp_fit.frequencies=frequencies
last_point= (harmonic_range[-1]*noramp_fit.nd_param.omega)+(noramp_fit.nd_param.omega*0.5)
plot_frequencies=frequencies[np.where(frequencies<last_point)]
noramp_fit.test_frequencies=plot_frequencies
noramp_fit.secret_data_time_series=current_results
#likelihood_func=noramp_fit.kaiser_filter(current_results)
noramp_fit.pass_extra_data(current_results, False)
noramp_fit.optim_list=['E_0','k_0', 'Ru','Cdl', 'gamma', 'omega']
param_boundaries=[[0, 1e-11],[1e-6, 1.0000001e-9]]
noramp_fit.define_boundaries(param_boundaries)
num_vals=20
gamma_vals=np.linspace(1e-12, 1e-9, num_vals)
cdl_vals=np.linspace(0, 1e-6, num_vals)
ru_vals=np.linspace(2.7, 1000, num_vals)
e0_vals=np.linspace(0.25, 0.6, num_vals)
k0_vals=np.linspace(0, 1e4, num_vals)
e0_vals=np.linspace(estart, ereverse, num_vals)
phase_vals=np.linspace(0, 2*math.pi, num_vals)
error_mat=np.zeros((num_vals, num_vals))
bounds_vals=np.arange(1, 20, 2)
print param_list['E_0']
means=[0.15880181408980204, 649.8257646670029, 1000.1079207626844, 0.057824979417721606e-6, 6.110952927721241e-10, 8.856800099582946]


#means=[2.85905314e-01, 5.86252081e+00, 1.19877032e-10, 4.19903721e-05, 4.78258908e+02, 8.94055732e+00]
#means2=[9.30117566e+03,1.16018837e-10,1.05469193e-05,1.63103403e+03,8.94085502e+00]
methods=["Bisect", "Brent minimisation", "Newton-Raphson"]
widths=[2,3,4]
alphas=[1, 0.8, 0.4]
ls=["-", "--", ":"]
for i in range(0, len(methods)):
    #for j in range(0, num_vals):
            #series=noramp_fit.simulate([cdl_vals[i], gamma_vals[j]], frequencies, "no_optimisation", "fourier")
            #start=np.where(noramp_fit.time_vec>0.012)
            #for j in range(int(0), len(noramp_fit.time_vec)):
            #for j in range(0, len(bounds_vals)):
                noramp_fit.bounds_val=100000
                noramp_fit.numerical_method=methods[i]
                #print bounds_vals[j]
                #noramp_fit.debug_time=noramp_fit.time_vec[j]
                start=time.time()
                time_series=noramp_fit.simulate(means, frequencies, "no_optimisation", "timeseries")
                end=time.time()-start
                plt.subplot(1,3,i+1)
                plt.plot(time_results, time_series)


                #plt.plot(time_results, time_series)
                #plt.plot(time_results, current_results)
                plt.title(methods[i]+" t="+str(end))

plt.show()

            #plt.plot(time_series)
            #plt.show()
            #time_series=np.array(time_series)
            #if len(time_series[np.where(abs(time_series)>200)])>0:
        #        print [cdl_vals[i], gamma_vals[j], ru_vals[k],k0_vals[lcv_4],phase_vals[lcv_5]],
    #            plt.plot(time_series)
    #            plt.show()
    #ax=sns.heatmap(error_mat, xticklabels=(gamma_vals), yticklabels=cdl_vals, fmt=".2g")
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
