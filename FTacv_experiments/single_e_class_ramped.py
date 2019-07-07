import isolver_ramped
import isolver_ramped_dispersion
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm
from params_class import params
import copy
import time
class single_electron:
    def __init__(self, dim_paramater_dictionary, simulation_options, other_values):
        key_list=dim_paramater_dictionary.keys()
        self.nd_param=params(dim_paramater_dictionary)
        #for i in range(0, len(key_list)):
        #    self.nd_param.non_dimensionalise(key_list[i], dim_paramater_dictionary[key_list[i]])
        self.nd_param=params(dim_paramater_dictionary)
        self.dim_dict=copy.deepcopy(dim_paramater_dictionary)
        self.nd_dict=self.nd_param.__dict__
        self.simulation_options=simulation_options
        self.optim_list=self.simulation_options["optim_list"]
        self.harmonic_range=other_values["harmonic_range"]
        self.num_harmonics=len(self.harmonic_range)
        self.filter_val=other_values["filter_val"]
        self.bounds_val=other_values["bounds_val"]
        if simulation_options["no_transient"]!=False:
            self.no_transient=simulation_options["no_transient"]
        if self.simulation_options["experimental_fitting"]==True:
            print "changing this"
            other_values["experiment_time"]=other_values["experiment_time"][:other_values["signal_length"]]/self.nd_param.c_T0
            other_values["experiment_current"]=other_values["experiment_current"][:other_values["signal_length"]]/self.nd_param.c_I0
            other_values["experiment_voltage"]=other_values["experiment_voltage"][:other_values["signal_length"]]/self.nd_param.c_E0
            self.time_vec=other_values["experiment_time"]
            if self.simulation_options["no_transient"]!=False:
                other_values["experiment_current"]=self.transient_remover(self.no_transient, other_values["experiment_time"], other_values["experiment_current"])
        else:
            self.nd_param.time_end=(self.nd_param.num_peaks/self.nd_param.nd_omega)*2*math.pi
            self.times(other_values["signal_length"])
        frequencies=np.fft.fftfreq(len(self.time_vec), self.time_vec[1]-self.time_vec[0])
        self.frequencies=frequencies[np.where(frequencies>0)]
        last_point= (self.harmonic_range[-1]*self.nd_param.omega)+(self.nd_param.omega*self.filter_val)
        self.test_frequencies=frequencies[np.where(self.frequencies<last_point)]
        self.other_values=other_values
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
    def kinetic_dispersion(self):
        k0_weights=np.zeros(self.simulation_options["dispersion_bins"])
        if (self.nd_param.k0_loc-(self.nd_param.k0_range/2))<0:
            k_start=0
            k_end=self.nd_param.k0_range
        else:
            k_start=self.nd_param.k0_loc-(self.nd_param.k0_range/2)
            k_end=self.nd_param.k0_loc+(self.nd_param.k0_range/2)
        k0_vals=np.linspace(k_start,k_end, self.simulation_options["dispersion_bins"])
        k0_weights[0]=lognorm.cdf(k0_vals[0], self.nd_param.k0_shape, loc=self.nd_param.k0_loc, scale=self.nd_param.k0_scale)
        for k in range(1, len(k0_weights)):
            k0_weights[k]=lognorm.cdf(k0_vals[k], self.nd_param.k0_shape, loc=self.nd_param.k0_loc, scale=self.nd_param.k0_scale)-lognorm.cdf(k0_vals[k-1], self.nd_param.k0_shape, loc=self.nd_param.k0_loc, scale=self.nd_param.k0_scale)
        #plt.plot(k0_vals, k0_weights)
        #plt.title("k0")
        #plt.show()
        return k0_vals, k0_weights
    def therm_dispersion(self):
        e0_weights=np.zeros(self.simulation_options["dispersion_bins"])
        e0_vals=np.linspace(self.nd_param.E_start,self.nd_param.E_reverse, self.simulation_options["dispersion_bins"])
        e0_weights[0]=norm.cdf(e0_vals[0], loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        for i in range(1, len(e0_weights)):
            e0_weights[i]=norm.cdf(e0_vals[i],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)-norm.cdf(e0_vals[i-1],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        #plt.plot(e0_vals, e0_weights)
        #print self.nd_param.E0_mean,self.nd_param.E0_std
        #plt.title("e0")
        #plt.show()
        return e0_vals, e0_weights
    def weight_matrix(self,e0_disp, k0_disp):
        e0_mat, k0_mat=np.meshgrid(e0_disp, k0_disp)
        weights=np.multiply(e0_mat, k0_mat)
        return weights
    def test_vals(self, parameters, likelihood, test=False):
        orig_likelihood=self.simulation_options["likelihood"]
        orig_label=self.simulation_options["label"]
        self.simulation_options["likelihood"]=likelihood
        self.simulation_options["label"]="MCMC"
        results=self.simulate(parameters, self.frequencies, test)
        self.simulation_options["likelihood"]=orig_likelihood
        self.simulation_options["label"]=orig_label
        return results
    def simulate(self,parameters, frequencies, test=False):
        if len(parameters)!= len(self.optim_list):
            raise ValueError('Wrong number of parameters')
        if self.simulation_options["label"]=="cmaes":
            normed_params=self.change_norm_group(parameters, "un_norm")
        else:
            normed_params=copy.deepcopy(parameters)
        for i in range(0, len(self.optim_list)):
            self.dim_dict[self.optim_list[i]]=normed_params[i]
            #print self.optim_list[i]
        self.nd_param=params(self.dim_dict)
        if self.simulation_options["dispersion"]==True:
            e0_vals, e0_disp=self.therm_dispersion()
            k0_vals, k0_disp=self.kinetic_dispersion()
            weights=self.weight_matrix(e0_disp, k0_disp)
            print np.sum(e0_disp), np.sum(k0_disp), np.sum(np.sum(weights))
            if np.sum(k0_disp)<0.99:
                print "FIRED"
                return np.zeros(len(self.time_vec))
            start=time.time()
            #time_series_disp=isolver_ramped_dispersion.e_surface_dispersed(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,e0_disp, k0_disp, weights, self.time_vec)
            disp_series=np.zeros(len(self.time_vec))
            time_series2=isolver_ramped.e_surface(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0, self.time_vec)
            for i in range(0, len(e0_vals)):
                for j in range(0, len(k0_vals)):
                    time_series=isolver_ramped.e_surface(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,e0_vals[i], k0_vals[j], self.time_vec)
                    disp_series=np.add(disp_series, np.multiply(time_series,weights[i][j]))
            plt.plot(time_series2)
            plt.plot(disp_series)
            plt.show()
            print time.time()-start, "TIME"
        else:
            time_series=isolver_ramped.e_surface(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0, self.time_vec)
        if self.simulation_options["no_transient"]==True:
            new_array=np.zeros(len(time_series))
            time_series=np.array(time_series)
            new_array[self.time_idx]=time_series[self.time_idx]
            time_series=new_array
        time_series=np.array(time_series)
        #time_series=np.flip(time_series)
        #time_series=np.flip(time_series, axis=0)
        #time_series=time_series*-1
        #time_series=np.flip(time_series*-1, axis=0)
        #time_series=(time_series*-1)
        if self.simulation_options["likelihood"]=='fourier':
            filtered=self.kaiser_filter(time_series)
            if (self.simulation_options["test"]==True or test==True):
                #plt.subplot(1,3,1)
                #plt.title("Likelihood function ("+str(self.harmonic_range[0])+"-"+str(self.harmonic_range[-1])+"harmonics)")
                plt.plot(self.secret_data_fourier, label="data")
                plt.plot(filtered , alpha=0.7, label="numerical")
                #plt.legend()
                #fourier_end=12
                #plt.title("Fourier spectrum up to harmonic " + str(fourier_end))
                #self.fourier_plotter(self.secret_data_time_series, "Data", time_series, "numerical", fourier_end)
                plt.show()

            return filtered
        elif self.simulation_options["likelihood"]=='timeseries':
            if self.simulation_options["test"]==True:

                plt.subplot(1,2,1)
                plt.plot(self.other_values["experiment_voltage"],time_series)
                plt.plot(self.other_values["experiment_voltage"],self.secret_data_time_series, alpha=0.7)
                plt.subplot(1,2,2)
                plt.plot(self.other_values["experiment_time"],time_series)
                plt.plot(self.other_values["experiment_time"],self.secret_data_time_series, alpha=0.7)
                plt.show()
            return (time_series)
