import isolver_noramp
import isolver_martin
import isolver_martin_newton
import isolver_martin_bisect
import isolver_brent
import isolver_inverted
from scipy.stats import norm, lognorm
import math
import numpy as np
import matplotlib.pyplot as plt
from params_class import params
from decimal import Decimal
import copy
import time
class single_electron:
    def __init__(self, dim_paramater_dictionary, simulation_options, other_values):
        key_list=list(dim_paramater_dictionary.keys())
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
            self.time_idx=simulation_options["no_transient"]
        else:
            self.time_idx=0
        if self.simulation_options["experimental_fitting"]==True:
            self.time_vec=other_values["experiment_time"][:other_values["signal_length"]]/self.nd_param.c_T0
            other_values["experiment_time"]=other_values["experiment_time"][self.time_idx:other_values["signal_length"]]/self.nd_param.c_T0
            other_values["experiment_current"]=other_values["experiment_current"][self.time_idx:other_values["signal_length"]]/self.nd_param.c_I0
            other_values["experiment_voltage"]=other_values["experiment_voltage"][self.time_idx:other_values["signal_length"]]/self.nd_param.c_E0

        else:
            self.nd_param.time_end=(self.nd_param.num_peaks/self.nd_param.nd_omega)*2*math.pi
            self.times(other_values["signal_length"])
        frequencies=np.fft.fftfreq(len(self.time_vec), self.time_vec[1]-self.time_vec[0])
        self.frequencies=frequencies[np.where(frequencies>0)]
        last_point= (self.harmonic_range[-1]*self.nd_param.omega)+(self.nd_param.omega*self.filter_val)
        self.test_frequencies=frequencies[np.where(self.frequencies<last_point)]
        self.other_values=other_values
        self.boundaries=None
        self.counter=0
    def define_boundaries(self, param_bounds):
        self.param_bounds=param_bounds
    def def_optim_list(self, optim_list):
        keys=list(self.dim_dict.keys())
        for i in range(0, len(optim_list)):
            if optim_list[i] in keys:
                continue
            else:
                raise KeyError("Parameter " + optim_list[i]+" not in list")
        self.optim_list=optim_list
        param_boundaries=np.zeros((2, self.n_parameters()))
        for i in range(0, self.n_parameters()):
                param_boundaries[0][i]=self.param_bounds[self.optim_list[i]][0]
                param_boundaries[1][i]=self.param_bounds[self.optim_list[i]][1]

        self.boundaries=param_boundaries

        if ("E0_std" in optim_list) or ("k0_loc" in optim_list):
            self.simulation_options["dispersion"]=True
        else:
            self.simulation_options["dispersion"]=False
        if "phase" in optim_list and "cap_phase" not in optim_list:
            self.phase_only=True
        else:
            self.phase_only=False
    def normalise(self, norm, boundaries):
        return  (norm-boundaries[0])/(boundaries[1]-boundaries[0])
    def un_normalise(self, norm, boundaries):
        return (norm*(boundaries[1]-boundaries[0]))+boundaries[0]
    def i_nondim(self, current):
        return np.multiply(current, self.nd_param.c_I0)
    def e_nondim(self, potential):
        return np.multiply(potential, self.nd_param.c_E0)
    def t_nondim(self, time):
        return np.multiply(time, self.nd_param.c_T0)
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
        #comp_results=np.append((np.real(results)), np.imag(results))
        return abs(results)
    def abs_transform(self, data):
        window=np.hanning(len(data))
        hanning_transform=np.multiply(window, data)
        f_trans=abs(np.fft.fft(hanning_transform[len(data)/2+1:]))
        return f_trans
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
                print(normed_params[i], self.optim_list[i],[self.boundaries[0][i],self.boundaries[1][i]])
                normed_params[i]=self.normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
        return normed_params
    def variable_returner(self):
        variables=vars(self.nd_param)
        for key in list(variables.keys()):
            if type(variables[key])==int or type(variables[key])==float or type(variables[key])==np.float64:
                print(key, variables[key])
    def pick_paramaters(self, param_vals, desired_params):
        if len(desired_params)>len(param_vals):
            raise ValueError("Too many parameters")
        num_params=len(desired_params)
        params=np.zeros(num_params)
        for i in range(0,num_params):
            idx=self.optim_list.index(desired_params[i])
            print(desired_params[i], idx)
            params[i]=param_vals[idx]
        return params
    def param_scanner(self, param_vals, param_list, unit_dict,percent, title, boundaries=False):

        unit_list=[unit_dict[k] for k in param_list]
        current_optim_list=self.optim_list
        self.optim_list=param_list
        num_params=len(param_list)
        for i in range(1, num_params):
            if num_params%i==0:
                col=i

        pc_change=[1-percent, 1, 1+percent]
        rows=num_params/col
        first_elem=((np.arange(0, rows)*col))
        bottom_elem=np.arange((rows*col)-col, rows*col)
        for i in range(0, num_params):
            ax=plt.subplot(rows, col, i+1)
            if i in first_elem:
                ax.set_ylabel("Current(A)")
            else:
                ax.axes.get_yaxis().set_ticks([])
            if i in bottom_elem:
                ax.set_xlabel("Voltage(V)")
            else:
                ax.axes.get_xaxis().set_ticks([])
            true_val=param_vals[i]
            plt.title(param_list[i])
            if boundaries==False:
                var_vals=np.multiply(true_val, pc_change)
            else:
                var_vals=[self.param_bounds[param_list[i]][0], true_val, self.param_bounds[param_list[i]][1]]
            for j in range(0,3):
                print(var_vals)
                param_vals[i]=var_vals[j]
                time_series=self.test_vals(param_vals, "timeseries", test=False)
                voltages=self.define_voltages()
                if abs(var_vals[j])<0.01:
                    label="{:.3E}".format(Decimal(str(var_vals[j])))
                else:
                    label="%.3f" % var_vals[j]
                plt.plot(self.e_nondim(self.other_values["experiment_voltage"]), self.i_nondim(time_series),label=label+" "+ str(unit_dict[param_list[i]]), alpha=0.7)#
            param_vals[i]=true_val
            plt.legend()
        plt.suptitle(title)
        self.optim_list=current_optim_list
        plt.show()

    def test_vals(self, parameters, likelihood, test=False):
        orig_likelihood=self.simulation_options["likelihood"]
        orig_label=self.simulation_options["label"]
        self.simulation_options["likelihood"]=likelihood
        self.simulation_options["label"]="MCMC"
        results=self.simulate(parameters, self.frequencies, test)
        self.simulation_options["likelihood"]=orig_likelihood
        self.simulation_options["label"]=orig_label
        return results
    def kinetic_dispersion(self):
        k0_weights=np.zeros(self.simulation_options["dispersion_bins"])
        k_start=0
        k_end=1e4
        k0_vals=np.linspace(k_start,k_end, self.simulation_options["dispersion_bins"])
        k0_weights[0]=lognorm.cdf(k0_vals[0], self.nd_param.k0_shape, loc=self.nd_param.k0_loc, scale=self.nd_param.k0_scale)
        for k in range(1, len(k0_weights)):
            k0_weights[k]=lognorm.cdf(k0_vals[k], self.nd_param.k0_shape, loc=self.nd_param.k0_loc, scale=self.nd_param.k0_scale)-lognorm.cdf(k0_vals[k-1], self.nd_param.k0_shape, loc=self.nd_param.k0_loc, scale=self.nd_param.k0_scale)
        #plt.plot(k0_vals, k0_weights)
        #plt.title("k0")
        #plt.show()
        return k0_vals, k0_weights

    def therm_dispersion(self):
        self.e0_min=self.nd_param.E_start#self.nd_param.E0_mean-(5*self.nd_param.E0_std)
        self.e0_max=self.nd_param.E_reverse#self.nd_param.E0_mean+(5*self.nd_param.E0_std)
        e0_weights=np.zeros(self.simulation_options["dispersion_bins"])
        e0_vals=np.linspace(self.e0_min, self.e0_max, self.simulation_options["dispersion_bins"])
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
    def numerical_plots(self):
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

    def weight_matrix(self,e0_disp, k0_disp):
        e0_mat, k0_mat=np.meshgrid(e0_disp, k0_disp)
        weights=np.multiply(e0_mat, k0_mat)
        return weights
    def variable_params(self, param_list, param_vals):
        true_vals=copy.deepcopy(param_list)
        y=len(param_list)#
        for i in range(1, y):
            if (y%i)==0:
                x=i
        for i in range(0, y):
            plt.subplot((y/x), x, i+1)


    def simulate(self,parameters, frequencies, test=False):
        if len(parameters)!= len(self.optim_list):
            print(self.optim_list)
            print(parameters)
            raise ValueError('Wrong number of parameters')
        if self.simulation_options["label"]=="cmaes":
            normed_params=self.change_norm_group(parameters, "un_norm")
        else:
            normed_params=copy.deepcopy(parameters)
        for i in range(0, len(self.optim_list)):
            self.dim_dict[self.optim_list[i]]=normed_params[i]
        if self.phase_only==True:
            self.dim_dict["cap_phase"]=self.dim_dict["phase"]
        self.nd_param=params(self.dim_dict)
        if self.simulation_options["numerical_method"]=="Bisect":
            solver=isolver_martin_bisect.martin_surface_bisect
        elif self.simulation_options["numerical_method"]=="Brent minimisation":
            solver=isolver_brent.martin_surface_brent
        elif self.simulation_options["numerical_method"]=="Newton-Raphson":
            solver=isolver_martin_newton.martin_surface_newton
        elif self.simulation_options["numerical_method"]=="inverted":
            solver=isolver_inverted.martin_surface_brent
        if self.simulation_options["numerical_debugging"]!=False:
            self.numerical_plots()
        else:
            if self.simulation_options["dispersion"]==True:
                time_series=np.zeros(len(self.time_vec))
                bins=self.simulation_options["dispersion_bins"]
                if ("E0_std" in self.optim_list) and ("k0_shape" in self.optim_list):
                    e0_vals, e0_disp=self.therm_dispersion()

                    k0_vals, k0_disp=self.kinetic_dispersion()
                    weights=self.weight_matrix(e0_disp, k0_disp)
                    print(sum(k0_disp))
                    for i in range(0, bins):
                        for j in range(0, bins):
                            time_series_current=solver(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,e0_vals[i], k0_vals[j],self.nd_param.cap_phase,self.time_vec[-1], self.time_vec, -1, self.bounds_val)
                            time_series=np.add(time_series, np.multiply(time_series_current, weights[i,j]))
                elif "E0_std" in self.optim_list:
                    e0_vals, e0_disp=self.therm_dispersion()
                    for i in range(0, bins):
                            time_series_current=solver(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,e0_vals[i], self.nd_param.k_0,self.nd_param.cap_phase,self.time_vec[-1], self.time_vec, -1, self.bounds_val)
                            time_series=np.add(time_series, np.multiply(time_series_current, e0_disp[i]))
                elif "k0_shape" in self.optim_list:
                    k0_vals, k0_disp=self.kinetic_dispersion()
                    for i in range(0, bins):
                        time_series_current=solver(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, k0_vals[i],self.nd_param.cap_phase,self.time_vec[-1], self.time_vec, -1, self.bounds_val)
                        time_series=np.add(time_series, np.multiply(time_series_current, k0_disp[i]))
                else:
                    time_series=solver(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0,self.nd_param.cap_phase,self.time_vec[-1], self.time_vec, -1, self.bounds_val)

            else:
                time_series=solver(self.nd_param.Cdl, self.nd_param.CdlE1, self.nd_param.CdlE2,self.nd_param.CdlE3, self.nd_param.nd_omega, self.nd_param.phase, math.pi,self.nd_param.alpha, self.nd_param.E_start,  self.nd_param.E_reverse, self.nd_param.d_E, self.nd_param.Ru, self.nd_param.gamma,self.nd_param.E_0, self.nd_param.k_0,self.nd_param.cap_phase, self.time_vec[-1], self.time_vec, -1, self.bounds_val)


        if self.simulation_options["no_transient"]!=True:
            time_series=time_series[self.time_idx:]
        time_series=np.array(time_series)
        if self.simulation_options["likelihood"]=='fourier':
            filtered=self.kaiser_filter(time_series)
            if (self.simulation_options["test"]==True or test==True):
                plt.plot(self.secret_data_fourier, label="data")
                plt.plot(filtered , alpha=0.7, label="numerical")
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
