from single_e_anal import analytical_electron
from single_e_class_anal_comparison import single_electron
import math
import numpy as np
import matplotlib.pyplot as plt
freq_range=[5*(10**x) for x in np.arange(0,6)]
time_ends=[10.0, 5.0,0.7, 0.07, 0.015, 0.015]
num_freqs=len(freq_range)
plt.figure(num=None, figsize=(9,10), dpi=120, facecolor='w', edgecolor='k')
for lcv_1 in range(0, num_freqs):
    param_list={
        'E_start': -0.5, #(starting dc voltage - V)
        'E_reverse': 0.5,    #  (reverse dc voltage - V)
        'omega':freq_range[lcv_1],#8.88480830076,  #    (frequency Hz)
        'd_E': 0.5,   #(ac voltage amplitude - V) freq_range[j],#
        'v': 150.97e-3,   #       (scan rate s^-1)
        'area': 0.03, #(electrode surface area cm^2)
        'Ru': 0.0,  #     (uncompensated resistance ohms)
        'Cdl': 1e-4*0, #(capacitance parameters)
        'CdlE1': 0,#0.000653657774506,
        'CdlE2': 0,#0.000245772700637,
        'CdlE3': 0,#1.10053945995e-06,
        'gamma': 6.5e-12,          # (surface coverage per unit area)
        'E_0': -0.1,      #       (reversible potential V)
        'E0_std':0.0312279186927,# (reversible potential dispersion)
        'k_0': 10.0, #(reaction rate s-1)
        'k0_std': 0.0,
        'alpha': 0.55,
        'sampling_freq' : (1.0/200),
        'phase' : 3*(math.pi/2)*0,
        'time_end':1000
    }
    anal=analytical_electron(param_list,0.0, 1.0)
    time_range=np.arange(0,time_ends[lcv_1], anal.nd_param.sampling_freq)
    non_dim_time_range=time_range*anal.nd_param.c_T0
    #non_dim_time_range=np.multiply(time_range, param_list['k_0'])
    time_len=len(time_range)
    i_val=np.zeros(time_len)-1
    for j in range(0, time_len-1):
        i_val[j]=anal.i(non_dim_time_range[j])

    nondim_i=anal.nd_param.c_I0
    #plt.subplot(2, 3, i+1)

    #f=np.fft.fftfreq(time_len, anal.nd_param.sampling_freq)
    #Y=np.fft.fft(i_val)
    harmonic_range=np.arange(1,10,1)
    numerical=single_electron(param_list, [], harmonic_range, 0.5)
    numerical.time_vec=non_dim_time_range
    numerical.num_points=time_len
    numerical.label="MCMC"
    synthetic_data=numerical.simulate([], time_range, "nah", "timeseries", "no")
    synthetic_data=np.multiply(synthetic_data, nondim_i)
    i_val=np.multiply(i_val, nondim_i)
    peak_ratio=(max(synthetic_data[time_len/2:]))/(max(i_val[time_len/2:]))
    peak_ratio_inv=1/peak_ratio
    #nd_dict=vars(param_list)
    #keys=nd_dict.keys()
    #for i in range(0, len(keys)):
    #    print keys[i], param_list[keys[i]]
        #print nd_dict[keys[i]]
    #plt.axhline(max(i_val[time_len/2:]), color="black", linestyle="--")
    #plt.axhline(max(synthetic_data[time_len/2:]), color="black", linestyle="--")
    plt.subplot(2, num_freqs/2, lcv_1+1)
    plt.plot(time_range, i_val*1000, label="analytical")#1.1245022593473897
    plt.plot(time_range, synthetic_data*1000, label="numerical")
    percent_diff=abs(peak_ratio-peak_ratio_inv)*100
    #plt.plot(np.subtract(i_val, synthetic_data))
    plt.title(str(freq_range[lcv_1])+ "    "+"$\Delta$=" + str(round(percent_diff,3))+ "%")
    plt.legend()
    plt.xlabel('Time(s)')
    plt.ylabel('Current(mA)')
plt.subplots_adjust(left=0.05,right=0.98, bottom=0.05, top=0.95, wspace=0.27)
plt.show()
