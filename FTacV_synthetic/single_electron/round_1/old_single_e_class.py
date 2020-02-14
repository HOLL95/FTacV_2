import isolver
import math
import numpy as np
import matplotlib.pyplot as plt
from params_class import params
import copy
import time
from scipy import signal
class single_electron:
    def __init__(self, dim_paramater_dictionary, optim_list, harmonic_range, filter_val):
        key_list=list(dim_paramater_dictionary.keys())
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
    def score_func(self, best_guess):
        score=sum(np.sqrt(np.power(np.subtract(best_guess, self.secret_data),2)))
        return score
    def kaiser_filter(self, time_series, harmonical=False):
        frequencies=self.frequencies
        L=len(time_series)
        window=np.hanning(L)
        time_series=np.multiply(time_series, window)
        f=np.fft.fftfreq(len(time_series), self.time_vec[1]-self.time_vec[0])
        Y=np.fft.fft(time_series)
        plt.plot(f,Y)
        plt.show()
        Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        top_hat=copy.deepcopy(Y[0:len(frequencies)])
        if harmonical==True:
            harmonics=np.zeros((self.num_harmonics, len(time_series)))
        for i in range(0, self.num_harmonics):
            true_harm=self.harmonic_range[i]*self.nd_param.omega*1.029
            #print true_harm, true_harm+(self.nd_param.omega*self.filter_val), true_harm-(self.nd_param.omega*self.filter_val)
            filter_bit=top_hat[np.where((frequencies<(true_harm+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm-(self.nd_param.omega*self.filter_val)))]
            #filter_bit=np.multiply(filter_bit, np.kaiser(len(filter_bit), 50))
            top_hat[np.where((frequencies<(true_harm+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm-(self.nd_param.omega*self.filter_val)))]=filter_bit
            if harmonical==True:
                print((self.harmonic_range[i]))
                harmonics[i,:][np.where((frequencies<(true_harm+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm-(self.nd_param.omega*self.filter_val)))]=filter_bit
                plt.plot((np.fft.ifft(harmonics[i,:])))
                plt.show()
        #if harmonical==True:
        #    return harmonics
        first_harm=self.harmonic_range[0]*self.nd_param.omega*1.029
        #top_hat=top_hat[np.where((frequencies>((self.harmonic_range[0]*self.nd_param.omega)-(self.nd_param.omega*0.5)))&(frequencies<((self.harmonic_range[-1]*self.nd_param.omega)+(self.nd_param.omega*0.5))))]
        plt.plot(top_hat)
        plt.show()
        #plt.plot(self.test_frequencies, top_hat)
        #plt.show()
        return (top_hat)
    def times(self, num_points):
        self.num_points=num_points
        self.time_vec=np.linspace(0, self.nd_param.time_end, num_points)
    def harmonics_plotter(self, original_harmonics, predicted_harmonics):
        fig, ax = plt.subplots(len(predicted_harmonics),1)
        for i in range(0, len(predicted_harmonics)):
            axes=ax[i]
            cuts=20000
            axes.plot(predicted_harmonics[i,:])#cuts:-cuts]
            axes.plot(original_harmonics[i,:], alpha=0.7)
        plt.show()
    def simulate(self,parameters, frequencies, flag='optimise', flag2='timeseries', test="no", score=False):

        normed_params=copy.deepcopy(parameters)
        if flag=='optimise' and self.label=="cmaes":
            for i in range(0,len(parameters)):
                normed_params[i]=self.un_normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
        for i in range(0, len(self.optim_list)):
                self.nd_param.non_dimensionalise(self.optim_list[i], normed_params[i])
        #results=isolver.I_tot_solver(self.nd_param.Cdl,self.nd_param.CdlE1,self.nd_param.CdlE2,self.nd_param.CdlE3,self.nd_param.nd_omega,self.nd_param.v, self.nd_param.alpha , \
        #                            self.nd_param.E_start,  self.nd_param.E_reverse,  self.nd_param.d_E,  self.nd_param.Ru, self.nd_param.sampling_freq, self.time_vec,  self.nd_param.gamma, \
        #                             self.nd_param.E_0, self.nd_param.k_0, self.nd_param.phase)
        time_series=isolver.e_surface(self.nd_param.Cdl,self.nd_param.CdlE1,self.nd_param.CdlE2,self.nd_param.CdlE3,self.nd_param.nd_omega,self.nd_param.v, self.nd_param.alpha , \
                                    self.nd_param.E_start,  self.nd_param.E_reverse,  self.nd_param.d_E,  self.nd_param.Ru,200, self.time_vec,  self.nd_param.gamma, \
                                     self.nd_param.E_0, self.nd_param.k_0, self.nd_param.phase, math.pi, self.num_points)

        #plt.plot(np.subtract(results[1],self.time_vec))
        #time_series=np.interp(self.time_vec,results[0], results[1])#desired co-oridnates, x-coords, y-coords
        if flag2=='fourier':
            filtered=self.kaiser_filter(time_series)
            if score ==True:
                return self.score_func(filtered)
            if test=="yes":
                print(normed_params)
                plt.plot(frequencies, self.secret_data_fourier)
                plt.plot(frequencies, filtered, alpha=0.7)
                plt.show()
                print(("score", score))
                #plt.plot(frequencies, filtered)
                #plt.plot(frequencies, self.secret_data)
                #plt.show()
                #original_harmonics=self.kaiser_filter(self.secret_data_time_series, frequencies, harmonical=True)
                #predicted_harmonics=self.kaiser_filter(time_series, frequencies, harmonical=True)
                #self.harmonics_plotter(original_harmonics, predicted_harmonics)

            return filtered
        elif flag2=='timeseries':
            if test=="yes":
                plt.plot(self.time_vec, time_series)
                plt.plot(self.time_vec, self.secret_data_time_series, alpha=0.7)
                plt.show()
            return time_series

"""

                    def top_hat_filter(self, time_series, frequencies, harmonical=False):
                        #next1=2**(math.ceil(math.log(len(output))/math.log(2)))
                        #pad=np.zeros(int(next1)-len(output))
                        #time_series=np.concatenate([pad,output])
                        L=len(time_series)
                        T=self.nd_param.sampling_freq
                        window=np.hanning(L)
                        time_series=np.multiply(time_series, window)
                        f=np.fft.fftfreq(len(time_series), self.time_vec[1]-self.time_vec[0])
                        Y=np.fft.fft(time_series)
                        Y_pow=np.power(Y[0:len(frequencies)],2)
                        top_hat=np.zeros((self.num_harmonics, len(frequencies)), dtype=complex)
                        if harmonical == True:
                            harmonics=np.zeros((self.num_harmonics, len(frequencies)), dtype=complex)
                        for i in range(0, self.num_harmonics):
                            true_harm=self.harmonic_range[i]*self.nd_param.omega*1.029
                            top_hat[i,:]=Y[0:len(frequencies)]
                            top_hat[i,:][np.where(frequencies>(true_harm+(self.nd_param.omega*self.filter_val)))]=0
                            top_hat[i,:][np.where(frequencies<(true_harm-(self.nd_param.omega*self.filter_val)))]=0
                            if harmonical==True:
                                harmonics[i,:]=abs(np.fft.ifft(top_hat[i,:]))
                        top_hat=abs((np.sum(top_hat, axis=0)))
                        plt.plot(f,Y)
                        plt.plot(frequencies, top_hat)
                        plt.show()
                        if harmonical==True:
                            return harmonics
                        return top_hat

"""
