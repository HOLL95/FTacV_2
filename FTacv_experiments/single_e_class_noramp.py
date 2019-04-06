import isolver_noramp
import isolver_martin
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
        plt.plot(freqs,abs(Y1), label=label1)
        plt.plot(freqs, abs(Y2), label=label2, alpha=0.7)
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
        #Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        top_hat=copy.deepcopy(Y[0:len(frequencies)])

        #for i in range(0, self.num_harmonics):
        #    true_harm=self.harmonic_range[i]*self.nd_param.omega*self.nd_param.c_T0
        #    filter_bit=top_hat[np.where((frequencies<(true_harm+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm-(self.nd_param.omega*self.filter_val)))]
        #    results[np.where((frequencies<(true_harm+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm-(self.nd_param.omega*self.filter_val)))]=filter_bit
        true_harm=self.nd_param.omega*self.nd_param.c_T0
        first_harm=(self.harmonic_range[0]*true_harm)-(self.nd_param.omega*self.filter_val)
        last_harm=(self.harmonic_range[-1]*true_harm)+(self.nd_param.omega*self.filter_val)
        likelihood=top_hat[np.where((frequencies>first_harm) & (frequencies<last_harm))]
        results=np.zeros(len(top_hat), dtype=complex)
        results[np.where((frequencies>first_harm) & (frequencies<last_harm))]=likelihood
        return abs(results)
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
                normed_params[i]=self.normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
        return normed_params
    def simulate(self,parameters, frequencies, flag='optimise', flag2='fourier', test="no"):
        if len(parameters)!= len(self.optim_list):
            raise ValueError('Wrong number of parameters')
        if flag=='optimise' and self.label=="cmaes":
            normed_params=self.change_norm_group(parameters, "un_norm")
        else:
            normed_params=copy.deepcopy(parameters)
        for i in range(0, len(self.optim_list)):
                self.nd_param.non_dimensionalise(self.optim_list[i], normed_params[i])
        time_series=isolver_noramp.e_surface(self.nd_param.Cdl,self.nd_param.CdlE1,self.nd_param.CdlE2,self.nd_param.CdlE3,self.nd_param.nd_omega,1, self.nd_param.alpha , \
                                    self.nd_param.E_start,  self.nd_param.E_reverse,  self.nd_param.d_E,  self.nd_param.Ru,200, self.time_vec,  self.nd_param.gamma, \
                                     self.nd_param.E_0, self.nd_param.k_0, self.nd_param.phase, math.pi, self.num_points)
        time_series_2=isolver_martin.martin_surface(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0,self.time_vec[-1], self.time_vec)
        print self.nd_param.Ru
        plt.plot(time_series_2)
        plt.title("Martin's method")
        plt.show()
        if ('method_label' in vars(self)):
            if ("time" in self.method_label):
                flag2="timeseries"
            if "cdl" in self.method_label:
                time_series=self.non_faradaic_filter(time_series)

        if flag2=='fourier':
            filtered=self.kaiser_filter(time_series)
            if test=="yes":
                plt.subplot(1,2,1)
                plt.title("Likelihood function (4-9th harmonics)")
                plt.plot(self.secret_data_fourier[:len(filtered)/15], label="data")
                plt.plot(filtered[:len(filtered)/15] , alpha=0.7, label="numerical")
                plt.legend()
                plt.subplot(1,2,2)
                fourier_end=12
                plt.title("Fourier spectrum up to harmonic " + str(fourier_end))
                self.fourier_plotter(self.secret_data_time_series, "Data", time_series, "numerical", fourier_end)
                plt.show()

            return filtered
        elif flag2=='timeseries':
            if test=="yes":
                print normed_params
                plt.plot(self.time_vec, time_series)
                plt.plot(self.time_vec, self.secret_data_time_series, alpha=0.7)
                plt.show()
            return time_series
