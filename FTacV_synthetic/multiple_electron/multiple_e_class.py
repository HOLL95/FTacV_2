import sys
import time
import matplotlib.pyplot as plt
from params_class import params
import math
import numpy as np
import seq_electron_noramp
import seq_electron_classical
import copy
import pints.plot
class multiple_electron:
    def __init__(self, param_dict, method, harmonic_range):
        self.method=method
        self.d_params=param_dict
        self.key_list=param_dict.keys()
        self.nd_param=params(param_dict)
        for i in range(0, len(self.key_list)):
            self.nd_param.non_dimensionalise(self.key_list[i], param_dict[self.key_list[i]])

        self.harmonic_range=harmonic_range
        self.num_harmonics=harmonic_range[-1]-harmonic_range[0]
    def time_params(self, time_end):
        self.time_vec=np.arange(0, time_end, self.nd_param.sampling_freq)
        return self.time_vec
    def optimisation_params(self, str_params):
        self.params_for_opt=str_params
    def n_parameters(self):
        print len(self.boundaries[0])
        return len(self.boundaries[0])-1
    def pass_extra_data(self, data):
        self.secret_data=data
    def n_outputs(self):
        return 1
    def define_boundaries(self, boundaries):
        self.boundaries=boundaries
    def normalise(self, norm, boundaries):
        return  (norm-boundaries[0])/(boundaries[1]-boundaries[0])
    def un_normalise(self, norm, boundaries):
        return (norm*(boundaries[1]-boundaries[0]))+boundaries[0]
    def top_hat_filter(self, frequencies, data):
        window=np.hanning(len(data))
        data=np.multiply(window, data)
        Y=np.fft.fft(data)
        Y=Y[0:len(frequencies)]
        Y=np.power(Y,2)
        top_hat=np.zeros((self.num_harmonics, len(frequencies)), dtype=complex)
        for i in range(0, self.num_harmonics):
            top_hat[i,:]=Y
            top_hat[i,:][np.where(frequencies>(self.harmonic_range[i]+0.05))]=0
            top_hat[i,:][np.where(frequencies<(self.harmonic_range[i]-0.05))]=0
        top_hat=abs(np.sum(top_hat, axis=0))
        return top_hat


    def simulate(self, parameters,  frequencies, gen_data=False, test=False):
        #print parameters
        if gen_data==False:
            normed_params=copy.deepcopy(parameters)
            for i in range(0,len(parameters)):
                normed_params[i]=self.un_normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
            for i in range((2*self.nd_param.num_species), len(self.params_for_opt)):
                self.nd_param.non_dimensionalise(self.params_for_opt[i], normed_params[i])
            self.nd_param.non_dimensionalise("E_0", normed_params[0:self.nd_param.num_species])
            self.nd_param.non_dimensionalise("k_0", normed_params[self.nd_param.num_species:2*self.nd_param.num_species])
        #print self.nd_param.k_0, self.nd_param.E_0
        if self.method =="noramp":
            #start=time.time()
            results=seq_electron_noramp.current_solver(self.nd_param.Cdl,self.nd_param.CdlE1,self.nd_param.CdlE2,self.nd_param.CdlE3,self.nd_param.nd_omega,self.nd_param.v, self.nd_param.alpha ,
                                            self.nd_param.E_start,  self.nd_param.E_reverse,  self.nd_param.d_E,  self.nd_param.Ru, self.nd_param.sampling_freq, len(self.time_vec),  self.nd_param.gamma, \
                                             self.nd_param.phase,self.nd_param.num_species, self.nd_param.k_0, self.nd_param.E_0)
            #print "simulat" ,time.time()-start
        elif self.method =="classical":
            results=seq_electron_classical.current_solver(self.nd_param.Cdl,self.nd_param.CdlE1,self.nd_param.CdlE2,self.nd_param.CdlE3,self.nd_param.nd_omega,self.nd_param.v, self.nd_param.alpha ,
                                            self.nd_param.E_start,  self.nd_param.E_reverse,  self.nd_param.d_E,  self.nd_param.Ru, self.nd_param.sampling_freq, len(self.time_vec),  self.nd_param.gamma, \
                                             self.nd_param.phase,self.nd_param.num_species, self.nd_param.k_0, self.nd_param.E_0)
        #start=time.time()
        time_series=np.interp(self.time_vec, results[1], results[0])
        #print "interp",time.time()-start
        if gen_data==True:
            return time_series
        #start= time.time()
        filtered=self.top_hat_filter(frequencies, time_series)
    #    print "filter",time.time()-start
        if test==True:
            #plt.plot(self.time_vec, time_series)

            plt.semilogy(frequencies, self.secret_data, alpha=0.5)
            plt.semilogy(frequencies, filtered)
            plt.show()
        return filtered
