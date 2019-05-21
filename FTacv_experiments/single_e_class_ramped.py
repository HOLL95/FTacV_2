import isolver_ramped
import math
import numpy as np
import matplotlib.pyplot as plt
from params_class import params
import copy
import time
class single_electron:
    def __init__(self, dim_paramater_dictionary, optim_list, harmonic_range, filter_val):
        key_list=dim_paramater_dictionary.keys()
        self.nd_param=params(dim_paramater_dictionary)
        for i in range(0, len(key_list)):
            self.nd_param.non_dimensionalise(key_list[i], dim_paramater_dictionary[key_list[i]])
        nd_dict=self.nd_param.__dict__
        self.optim_list=optim_list
        self.harmonic_range=harmonic_range
        self.num_harmonics=len(harmonic_range)
        self.filter_val=filter_val
    def define_boundaries(self, boundaries):
        self.boundaries=boundaries
    def normalise(self, norm, boundaries):
        return  (norm-boundaries[0])/(boundaries[1]-boundaries[0])
    def un_normalise(self, norm, boundaries):
        return (norm*(boundaries[1]-boundaries[0]))+boundaries[0]
    def dimensionalise_data(self,data, times):
        data=np.multiply(data, self.nd_param.c_I0)
        times=np.multiply(times,self.nd_param.c_T0)
        return times, data
    def nondim_data(self, data, times):
        data=np.divide(data, self.nd_param.c_I0)
        times=np.divide(times,self.nd_param.c_T0)
        return dtimes, data
    def n_outputs(self):
        return 1
    def n_parameters(self):
        return len(self.optim_list)
    def pass_extra_data(self, time_series, fourier):
        self.secret_data_time_series=time_series
        self.secret_data_fourier=fourier
    def non_faradaic_filter(self, time_series):
        signal_length=len(time_series)
        filtered_val=np.zeros(signal_length)
        divisor=11
        endl=divisor-1
        filtered_val[:int(signal_length/divisor)]=time_series[:int(signal_length/divisor)]
        filtered_val[int(endl*signal_length/divisor):]=time_series[int(endl*signal_length/divisor):]
        return filtered_val
    def fourier_plotter(self, time_series_1, label1, time_series_2, label2, fourier_end):
        L=len(time_series_1)
        window=np.hanning(L)
        time_series_1=np.multiply(time_series_1, window)
        time_series_2=np.multiply(time_series_2, window)
        f=np.fft.fftfreq(len(time_series_1), self.time_vec[1]-self.time_vec[0])
        true_harm=self.nd_param.omega*self.nd_param.c_T0
        last_harm=(fourier_end*true_harm)+(self.nd_param.omega*self.filter_val)
        freqs=f[np.where((f>0) & (f<last_harm))]
        Y1=np.fft.fft(time_series_1)
        Y2=np.fft.fft(time_series_2)
        Y1=Y1[np.where((f>0) & (f<last_harm))]
        Y2=Y2[np.where((f>0) & (f<last_harm))]
        plt.subplot(1,3,2)
        plt.title("Real")
        plt.plot(freqs,np.real(Y1), label=label1)
        plt.plot(freqs, np.real(Y2), label=label2, alpha=0.7)
        plt.subplot(1,3,3)
        plt.title("Imaginary")
        plt.plot(freqs,np.imag(Y1), label=label1)
        plt.plot(freqs, np.imag(Y2), label=label2, alpha=0.7)
        plt.xlabel("frequency")
        plt.legend()
    def kaiser_filter(self, time_series, harmonical=False):
        frequencies=self.frequencies
        L=len(time_series)
        window=np.hanning(L)
        time_series=np.multiply(time_series, window)
        f=np.fft.fftfreq(len(time_series), self.time_vec[1]-self.time_vec[0])
        Y=np.fft.fft(time_series)
        #plt.plot(f,Y)
        #plt.show()
        #Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        top_hat=copy.deepcopy(Y[0:len(frequencies)])
        true_harm=self.nd_param.omega*self.nd_param.c_T0
        first_harm=(self.harmonic_range[0]*true_harm)-(self.nd_param.omega*self.filter_val)
        last_harm=(self.harmonic_range[-1]*true_harm)+(self.nd_param.omega*self.filter_val)
        likelihood=top_hat[np.where((frequencies>first_harm) & (frequencies<last_harm))]
        results=np.zeros(len(top_hat), dtype=complex)
        results[np.where((frequencies>first_harm) & (frequencies<last_harm))]=likelihood
        #comp_results=np.append(np.real(results), np.imag(results))
        return (results)
    def times(self, num_points):
        self.num_points=num_points
        #self.time_vec=np.arange(0, self.nd_param.time_end, self.nd_param.sampling_freq)
        self.time_vec=np.linspace(0, self.nd_param.time_end, num_points)
    def change_norm_group(self, param_list, method):
        normed_params=copy.deepcopy(param_list)
        if method=="un_norm":
            for i in range(0,len(param_list)):
                normed_params[i]=self.un_normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
        elif method=="norm":
            for i in range(0,len(param_list)):
                print normed_params[i], self.optim_list[i],[self.boundaries[0][i],self.boundaries[1][i]]
                normed_params[i]=self.normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
        return normed_params
    def variable_returner(self):
        variables=vars(self.nd_param)
        for key in variables.keys():
            if type(variables[key])==int or type(variables[key])==float or type(variables[key])==np.float64:
                print key, variables[key], type(variables[key])
    def transient_remover(self, start_time, times, current):
        time_idx=np.where(times>start_time)
        self.time_idx=time_idx
        new_array=np.zeros(len(current))
        new_array[self.time_idx]=current[self.time_idx]
        return  new_array

    def simulate(self,parameters, frequencies, flag='optimise', flag2='fourier', test="no"):
        var_list=vars(self)
        if len(parameters)!= len(self.optim_list):
            raise ValueError('Wrong number of parameters')
        if flag=='optimise' and self.label=="cmaes":
            normed_params=self.change_norm_group(parameters, "un_norm")
        else:
            normed_params=copy.deepcopy(parameters)
        #print normed_params
        for i in range(0, len(self.optim_list)):
                self.nd_param.non_dimensionalise(self.optim_list[i], normed_params[i])
        time_series=isolver_ramped.e_surface(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0, self.time_vec)
        if "no_transient" in var_list:
            new_array=np.zeros(len(time_series))
            time_series=np.array(time_series)
            new_array[self.time_idx]=time_series[self.time_idx]
            time_series=new_array
        time_series=np.flip(time_series)
        if flag2=='fourier':
            filtered=self.kaiser_filter(time_series)
            if test=="yes":
                #plt.subplot(1,3,1)
                plt.title("Likelihood function ("+str(self.harmonic_range[0])+"-"+str(self.harmonic_range[-1])+"harmonics)")
                plt.plot(self.secret_data_fourier)
                plt.plot(filtered, alpha=0.7)
                plt.legend()
                plt.show()
                #fourier_end=12
                #plt.title("Fourier spectrum up to harmonic " + str(fourier_end))
                #self.fourier_plotter(self.secret_data_time_series, "Data", time_series, "numerical", fourier_end)
                #plt.show()

            return filtered
        elif flag2=='timeseries':
            if test=="yes":
                plt.plot(time_series)
                plt.plot(self.secret_data_time_series, alpha=0.7)
                plt.show()
            return time_series
