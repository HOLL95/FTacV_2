import isolver_noramp
import isolver_martin
import isolver_martin_newton
import isolver_martin_bisect
import isolver_martin_brent
import isolver_inverted
import math
import numpy as np
import matplotlib.pyplot as plt
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
        print "OMEGA", self.nd_param.CdlE1
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
    def define_voltages(self):
        voltages=np.zeros(len(self.time_vec))
        for i in range(0, len(self.time_vec)):
            voltages[i]=isolver_noramp.et(self.nd_param.E_start, self.nd_param.E_reverse,self.nd_param.nd_omega, self.nd_param.phase, (self.time_vec[i]))
        return voltages
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
        time_series_2=np.multiply(time,e_series_2, window)
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
        #Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        top_hat=copy.deepcopy(Y[0:len(frequencies)])

        true_harm=self.nd_param.omega*self.nd_param.c_T0
        if sum(np.diff(self.harmonic_range))!=len(self.harmonic_range)-1:
            results=np.zeros(len(top_hat), dtype=complex)
            for i in range(0, self.num_harmonics):
                true_harm_n=true_harm*self.harmonic_range[i]
                index=[np.where((frequencies<(true_harm_n+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm_n-(self.nd_param.omega*self.filter_val)))]
                filter_bit=top_hat[index]
                results[index]=filter_bit
        else:
            first_harm=(self.harmonic_range[0]*true_harm)-(self.nd_param.omega*self.filter_val)
            last_harm=(self.harmonic_range[-1]*true_harm)+(self.nd_param.omega*self.filter_val)
            likelihood=top_hat[np.where((frequencies>first_harm) & (frequencies<last_harm))]
            results=np.zeros(len(top_hat), dtype=complex)
            results[np.where((frequencies>first_harm) & (frequencies<last_harm))]=likelihood
        comp_results=np.append(np.real(results), np.imag(results))
        return (comp_results)
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
        time_end=self.time_vec[-1]-start_time
        time_idx=np.where((times>start_time) & (times<time_end))
        self.time_idx=time_idx
        new_array=np.zeros(len(current))
        new_array[self.time_idx]=current[self.time_idx]
        return  new_array

    def simulate(self,parameters, frequencies, normalise =False, likelihood=None, test=False):
        if len(parameters)!= len(self.optim_list):
            raise ValueError('Wrong number of parameters')
        if self.simulation_options["label"]=="cmaes" or normalise==True:
            normed_params=self.change_norm_group(parameters, "un_norm")
        else:
            normed_params=copy.deepcopy(parameters)
        for i in range(0, len(self.optim_list)):
            self.dim_dict[self.optim_list[i]]=normed_params[i]
        self.nd_param=params(self.dim_dict)
        if self.simulation_options["numerical_method"]=="Bisect":
            solver=isolver_martin_bisect.martin_surface_bisect
        elif self.simulation_options["numerical_method"]=="Brent minimisation":
            solver=isolver_martin_brent.martin_surface_brent
        elif self.simulation_options["numerical_method"]=="Newton-Raphson":
            solver=isolver_martin_newton.martin_surface_newton
        elif self.simulation_options["numerical_method"]=="inverted":
            solver=isolver_inverted.martin_surface_brent
        if self.simulation_options["numerical_debugging"]!=False:
            self.debug_time=self.simulation_options["numerical_debugging"]
            time_series=solver(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0,self.time_vec[-1], self.time_vec, self.debug_time, self.bounds_val)
            current=time_series[0]
            residual=time_series[1]
            residual_gradient=time_series[2]
            #plt.subplot(1,2,1)
            #plt.semilogy(current, np.abs(residual))
            Nt=self.time_vec[-1]/self.nd_param.sampling_freq
            bounds_val=max(10*self.nd_param.Cdl*self.nd_param.d_E*self.nd_param.nd_omega/Nt,1.0)
            middle_index=(len(time_series[0])-1)/2 + 1
            I0=residual[middle_index]
            if  self.simulation_options["numerical_method"]=="Newton-Raphson":
                plt.subplot(1,2,1)
                plt.title("Residual, t="+str(self.debug_time))
                plt.plot(current, residual)
                plt.axvline(time_series[3][1], color="red",linestyle="--")
                plt.axvline(time_series[3][0]+time_series[3][2], color="black", linestyle="--")
                plt.axvline(time_series[3][0]-time_series[3][2], color="black",linestyle="--")
                plt.subplot(1,2,2)
                plt.title("Residual gradient")
                plt.plot(current, ((residual_gradient)))
                plt.show()
            else:
                plt.plot(current, residual, label=round(self.debug_time, 4))
                plt.show()
                #plt.axvline(time_series[3][0]+time_series[3][2], color="black",linestyle="--")
                #plt.axvline(time_series[3][0]-time_series[3][2],color="black",linestyle="--")

        else:

            time_series=solver(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0,self.time_vec[-1], self.time_vec, -1, self.bounds_val)
        if self.simulation_options["no_transient"]==True:
            new_array=np.zeros(len(time_series))
            time_series=np.array(time_series)
            new_array[self.time_idx]=time_series[self.time_idx]
            time_series=new_array
        time_series=np.array(time_series)
        time_series=np.flip(time_series)
        #time_series=np.flip(time_series, axis=0)
        #time_series=time_series*-1
        #time_series=np.flip(time_series*-1)
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
                plt.plot(time_series)
                plt.plot(self.secret_data_time_series, alpha=0.7)
                plt.show()
            return (time_series)
