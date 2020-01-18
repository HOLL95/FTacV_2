import isolver_martin_brent
#import isolver_martin_NR
from scipy.stats import norm, lognorm
import math
import numpy as np
import itertools
import multiprocessing as mp
import matplotlib.pyplot as plt
from params_class import params
from decimal import Decimal
import copy
import time
import pickle
import warnings
class single_electron:
    def __init__(self,file_name="", dim_paramater_dictionary={}, simulation_options={}, other_values={}, param_bounds={}, results_flag=True):
        if type(file_name) is dict:
            raise TypeError("Need to define a filename - this is currently a dictionary!")

        if len(dim_paramater_dictionary)==0 and len(simulation_options)==0 and len(other_values)==0:
            self.file_init=True
            file=open(file_name, "rb")
            save_dict=pickle.load(file, encoding="latin1")
            dim_paramater_dictionary=save_dict["param_dict"]
            simulation_options=save_dict["simulation_opts"]
            other_values=save_dict["other_vals"]
            param_bounds=save_dict["bounds"]
            self.save_dict=save_dict
        else:
            self.file_init=False

        key_list=list(dim_paramater_dictionary.keys())
        #for i in range(0, len(key_list)):
        #    self.nd_param.non_dimensionalise(key_list[i], dim_paramater_dictionary[key_list[i]])
        if simulation_options["method"]=="ramped":
            dim_paramater_dictionary["v_nondim"]=True
        self.nd_param=params(dim_paramater_dictionary)
        self.nd_param_dict=self.nd_param.nd_param_dict
        self.dim_dict=copy.deepcopy(dim_paramater_dictionary)
        self.nd_dict=self.nd_param.__dict__
        self.simulation_options=simulation_options
        self.optim_list=self.simulation_options["optim_list"]
        self.harmonic_range=other_values["harmonic_range"]
        self.num_harmonics=len(self.harmonic_range)
        self.filter_val=other_values["filter_val"]
        self.bounds_val=other_values["bounds_val"]
        self.time_array=[]
        if self.simulation_options["experimental_fitting"]==True:
            if simulation_options["method"]=="sinusoidal":
                time_end=(self.nd_param.num_peaks/self.nd_param.omega)
            elif simulation_options["method"]=="ramped":
                time_end=2*(self.nd_param.E_reverse-self.nd_param.E_start)
            elif simulation_options["method"]=="dcv":
                time_end=2*(self.nd_param.E_reverse-self.nd_param.E_start)
            if simulation_options["no_transient"]!=False:
                if simulation_options["no_transient"]>time_end:
                    warnings.warn("Previous transient removal method detected")
                    time_idx=tuple(np.where(other_values["experiment_time"]<time_end))
                    desired_idx=tuple((range(simulation_options["no_transient"],time_idx[0][-1])))
                    self.time_idx=desired_idx[0]
                else:
                    desired_idx=tuple(np.where((other_values["experiment_time"]<time_end) & (other_values["experiment_time"]>simulation_options["no_transient"])))
                    time_idx=tuple(np.where(other_values["experiment_time"]<time_end))
                    self.time_idx=desired_idx[0][0]
            else:
                desired_idx=tuple(np.where(other_values["experiment_time"]<time_end))
                time_idx=desired_idx
                self.time_idx=0
            if self.file_init==False or results_flag==True:
                self.time_vec=other_values["experiment_time"][time_idx]/self.nd_param.c_T0
                other_values["experiment_time"]=other_values["experiment_time"][desired_idx]/self.nd_param.c_T0
                other_values["experiment_current"]=other_values["experiment_current"][desired_idx]/self.nd_param.c_I0
                other_values["experiment_voltage"]=other_values["experiment_voltage"][desired_idx]/self.nd_param.c_E0
            else:
                if simulation_options["method"]=="sinusoidal":
                    self.nd_param.time_end=(self.nd_param.num_peaks)#/self.nd_param.omega)
                else:
                    self.nd_param.time_end=2*(self.nd_param.E_reverse-self.nd_param.E_start)
                self.times()

            #else:
            #    sf=other_values["experiment_time"][1]-other_values["experiment_time"][0]
            #    self.time_vec=np.arange(0, other_values["experiment_time"][-1], sf)



        else:
            if simulation_options["method"]=="sinusoidal":
                self.nd_param.time_end=(self.nd_param.num_peaks)#/self.nd_param.omega)
            else:
                self.nd_param.time_end=2*(self.nd_param.E_reverse-self.nd_param.E_start)
            self.times()
            if simulation_options["no_transient"]!=False:
                    transient_time=self.t_nondim(self.time_vec)
                    start_idx=np.where(transient_time>simulation_options["no_transient"])
                    self.time_idx=start_idx[0][0]
            else:
                    self.time_idx=0
        frequencies=np.fft.fftfreq(len(self.time_vec), self.time_vec[1]-self.time_vec[0])
        self.frequencies=frequencies[np.where(frequencies>0)]
        last_point= (self.harmonic_range[-1]*self.nd_param.omega)+(self.nd_param.omega*self.filter_val)
        self.test_frequencies=frequencies[np.where(self.frequencies<last_point)]
        self.other_values=other_values
        self.boundaries=None
        self.param_bounds=param_bounds
        if self.simulation_options["experimental_fitting"]==True:
            self.secret_data_fourier=self.kaiser_filter(other_values["experiment_current"])
            self.secret_data_time_series=other_values["experiment_current"]
    def define_boundaries(self, param_bounds):
        self.param_bounds=param_bounds
    def def_optim_list(self, optim_list):
        keys=list(self.dim_dict.keys())
        for i in range(0, len(optim_list)):
            if optim_list[i] in keys:
                continue
            else:
                raise KeyError("Parameter " + optim_list[i]+" not found in model")
        self.optim_list=optim_list
        param_boundaries=np.zeros((2, self.n_parameters()))
        check_for_bounds=vars(self)
        if "param_bounds" in list(check_for_bounds.keys()):
            for i in range(0, self.n_parameters()):
                    param_boundaries[0][i]=self.param_bounds[self.optim_list[i]][0]
                    param_boundaries[1][i]=self.param_bounds[self.optim_list[i]][1]

            self.boundaries=param_boundaries

        if ("E0_std" in optim_list) or ("k0_shape" in optim_list):
            self.simulation_options["dispersion"]=True
        else:
            self.simulation_options["dispersion"]=False
        if self.simulation_options["dispersion"]==True:
            if "k_0" not in self.optim_list:
                k0params=set(["k0_shape", "k0_scale"])
                if k0params.issubset(set(self.optim_list))==False:
                    missing_param=k0params-(k0params & set(self.optim_list))
                    raise ValueError("Missing the following kinetic dispersion parameters: "+ (", ").join([str(x) for x in missing_param]))
            if "E_0" not in self.optim_list:
                e0params=set(["E0_std", "E0_mean"])
                if e0params.issubset(set(self.optim_list))==False:
                    missing_param=e0params-(e0params & set(self.optim_list))
                    raise ValueError("Missing the following thermodynamic dispersion parameters: " + (", ").join([str(x) for x in missing_param]))
        if "phase" in optim_list and "cap_phase" not in optim_list:
            self.simulation_options["phase_only"]=True
        else:
            self.simulation_options["phase_only"]=False
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
    def define_voltages(self, transient=False):
        voltages=np.zeros(len(self.time_vec))
        if self.simulation_options["method"]=="sinusoidal":
            for i in range(0, len(self.time_vec)):
                voltages[i]=isolver_martin_brent.et(self.nd_param.E_start,self.nd_param.nd_omega, self.nd_param.phase, self.nd_param.d_E, (self.time_vec[i]))
        elif self.simulation_options["method"]=="ramped":
            for i in range(0, len(self.time_vec)):
                voltages[i]=isolver_martin_brent.c_et(self.nd_param.E_start, self.nd_param.E_reverse, (self.nd_param.E_reverse-self.nd_param.E_start) ,self.nd_param.nd_omega, self.nd_param.phase, 1,self.nd_param.d_E,(self.time_vec[i]))
        elif self.simulation_options["method"]=="dcv":
            for i in range(0, len(self.time_vec)):
                voltages[i]=isolver_martin_brent.dcv_et(self.nd_param.E_start, self.nd_param.E_reverse, (self.nd_param.E_reverse-self.nd_param.E_start) , 1,(self.time_vec[i]))
        #print(self.nd_param.E_start, self.nd_param.E_reverse, (self.nd_param.E_reverse-self.nd_param.E_start) ,self.nd_param.nd_omega, self.nd_param.phase, 1,self.nd_param.d_E,(self.time_vec[i]))
        if transient==True:
            voltages=voltages[self.time_idx:]
        return voltages
    def pass_extra_data(self, time_series, fourier):
        self.secret_data_time_series=time_series
        self.secret_data_fourier=fourier
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
            #self.test_frequencies=frequencies[np.where((frequencies>first_harm) & (frequencies<last_harm))]
            results=np.zeros(len(top_hat), dtype=complex)
            results[np.where((frequencies>first_harm) & (frequencies<last_harm))]=likelihood
        comp_results=np.append((np.real(results)), np.imag(results))
        return comp_results
    def abs_transform(self, data):
        window=np.hanning(len(data))
        hanning_transform=np.multiply(window, data)
        f_trans=abs(np.fft.fft(hanning_transform[len(data)/2+1:]))
        return f_trans
    def saved_param_simulate(self, params):
        if self.file_init==False:
            raise ValueError('No file provided')
        else:
            self.def_optim_list(self.save_dict["optim_list"])
            type=self.simulation_options["likelihood"]
            return self.test_vals(params,type, test=False)
    def save_state(self, results, filepath, filename, params):
        other_vals_save=self.other_values
        other_vals_save["experiment_time"]=results["experiment_time"]
        other_vals_save["experiment_current"]=results["experiment_current"]
        other_vals_save["experiment_voltage"]=results["experiment_voltage"]
        file=open(filepath+"/"+filename, "wb")
        save_dict={"simulation_opts":self.simulation_options, \
                    "other_vals":other_vals_save, \
                    "bounds":self.param_bounds, \
                    "param_dict":self.dim_dict ,\
                    "params":params, "optim_list":self.optim_list}
        pickle.dump(save_dict, file, pickle.HIGHEST_PROTOCOL)
        file.close()
    def calc_theta(self, current):
        voltages=self.define_voltages()
        if self.simulation_options["no_transient"]!=True:
            voltages=voltages[self.time_idx:]
        theta=np.zeros(len(current))
        theta[0]=0
        dt=self.nd_param.sampling_freq
        for i in range(1, len(current)):
            Er=voltages[i]-self.nd_param.Ru*current[i]
            expval1=Er-self.nd_param.E_0
            exp11=np.exp((1-self.nd_param.alpha)*expval1)
            exp12=np.exp((-self.nd_param.alpha)*expval1)
            u1n1_top=dt*self.nd_param.k_0*exp11 + theta[i-1]
            denom = ((dt*self.nd_param.k_0*exp11) +(dt*self.nd_param.k_0*exp12) + 1)
            theta[i]=u1n1_top/denom
        return theta
    def times(self):
        self.time_vec=np.arange(0, self.nd_param.time_end, self.nd_param.sampling_freq)
        #self.time_vec=np.linspace(0, self.nd_param.time_end, num_points)
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
    def param_scanner(self, param_list, unit_dict, title, num_scans, param_vals=[], pc=0):
        boundaries=self.param_bounds
        if len(param_vals)==0 and pc==0:
            pc_shift=False
            param_vals=[(boundaries[x][0]+boundaries[x][1])/2.0 for x in param_list]
        else:
            pc_shift=True
        unit_list=[unit_dict[k] for k in param_list]
        current_optim_list=self.optim_list
        self.optim_list=param_list
        num_params=len(param_list)
        for i in range(1, num_params):
            if num_params%i==0:
                col=i
        pc_change=np.arange(1-(0.05*num_scans/2), 1+(0.05*num_scans/2), 0.05)
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
            print(param_list[i])
            if pc_shift==True:
                var_vals=np.multiply(true_val, pc_change)
            else:
                var_vals=np.linspace(self.param_bounds[param_list[i]][0], self.param_bounds[param_list[i]][1], num_scans)
            for j in range(0,num_scans):
                param_vals[i]=var_vals[j]
                time_series=self.test_vals(param_vals, "timeseries", test=False)
                if type(self.other_values["experiment_voltage"])!=bool:
                    voltages=self.other_values["experiment_voltage"]
                else:
                    voltages=self.define_voltages()
                if abs(var_vals[j])<0.01:
                    label="{:.3E}".format(Decimal(str(var_vals[j])))
                else:
                    label="%.3f" % var_vals[j]
                plt.plot(self.e_nondim(voltages), self.i_nondim(time_series),label=label+" "+ str(unit_dict[param_list[i]]), alpha=0.7)#
            param_vals[i]=true_val
            plt.legend()
        plt.suptitle(title)
        self.optim_list=current_optim_list
        plt.show()

    def test_vals(self, parameters, likelihood, test=False):
        orig_likelihood=self.simulation_options["likelihood"]
        orig_label=self.simulation_options["label"]
        orig_test=self.simulation_options["test"]
        self.simulation_options["likelihood"]=likelihood
        self.simulation_options["label"]="MCMC"
        self.simulation_options["test"]=test
        results=self.simulate(parameters, self.frequencies)
        self.simulation_options["likelihood"]=orig_likelihood
        self.simulation_options["label"]=orig_label
        self.simulation_options["test"]=orig_test
        return results
    def paralell_disperse(self, solver):
        time_series=np.zeros(len(self.time_vec))
        if ("E0_std" in self.optim_list) and ("k0_shape" in self.optim_list):
            e0_vals, e0_disp=self.therm_dispersion()
            k0_vals, k0_disp=self.kinetic_dispersion()
            values=list(itertools.product(e0_vals, k0_vals))
            flags=list(zip(["E_0"]*len(values), ["k_0"]*len(values)))
            weights=list(itertools.product(e0_disp, k0_disp))
            weights=[weights[i][0]*weights[i][1] for i in range(len(weights))]
            weight_val_tuple=list(zip(flags, values, weights))
            paralell=paralell_class(self.nd_param_dict, self.time_vec, "sinusoidal", self.bounds_val, isolver_martin_brent.brent_current_solver)
            time_series=paralell.paralell_dispersion(weight_val_tuple)

        elif "E0_std" in self.optim_list:
            start1=time.time()
            e0_vals, e0_disp=self.therm_dispersion()

            if "alpha_dispersion" in self.simulation_options:
                alpha_bins=5
                alpha_vals=np.linspace(self.param_bounds["alpha"][0],self.param_bounds["alpha"][1], alpha_bins)
                alpha_weights=[1/alpha_bins]*alpha_bins
                for i in range(0, alpha_bins):
                    self.nd_param_dict["alpha"]=float(alpha_vals[i])
                    for j in range(0,self.simulation_options["dispersion_bins"]):
                        self.nd_param_dict["E_0"]=float(e0_vals[j])
                        time_series_current=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], -1, self.bounds_val)
                        time_series=np.add(time_series, np.multiply(time_series_current, e0_disp[j]*alpha_weights[i]))
            else:
                for i in range(0, self.simulation_options["dispersion_bins"]):
                    self.nd_param_dict["E_0"]=float(e0_vals[i])
                    time_series_current=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], -1, self.bounds_val)
                    time_series=np.add(time_series, np.multiply(time_series_current, e0_disp[i]))
        elif ("k0_shape" in self.optim_list):
            k0_vals, k0_disp=self.kinetic_dispersion()
            for i in range(0, self.simulation_options["dispersion_bins"]):
                self.nd_param_dict["k_0"]=k0_vals[i]
                time_series_current=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], -1, self.bounds_val)
                time_series=np.add(time_series, np.multiply(time_series_current, k0_disp[i]))
        print(self.dim_dict["k0_shape"], self.dim_dict["k0_scale"])
        #plt.plot(k0_vals, k0_disp)
        #plt.show()


        return time_series

    def generic_dispersion(self, solver):
        if "dispersion_parameters" not in self.simulation_options:
            raise ValueError("Dispersion parameters not defined")
        if len(self.simulation_options["dispersion_bins"])!=len(self.simulation_options["dispersion_parameters"]):
            raise ValueError("Need to define number of bins for each parameter")
        if len(self.simulation_options["dispersion_distributions"])!=len(self.simulation_options["dispersion_parameters"]):
            raise ValueError("Need to define distributions for each parameter")
        weight_arrays=[]
        weight_values=[]
        nd_dict=vars(self.nd_param)
        for i in range(0, len(self.simulation_options["dispersion_parameters"])):
            if self.simulation_options["dispersion_distributions"][i]=="uniform":
                if (self.simulation_options["dispersion_parameters"][i]+"_lower" not in self.dim_dict) or (self.simulation_options["dispersion_parameters"][i]+"_upper" not in self.dim_dict):
                    raise ValueError("Uniform distribution requires "+self.simulation_options["dispersion_parameters"][i]+"_lower and " + self.simulation_options["dispersion_parameters"][i]+"_upper")
                else:
                    weight_values.append(np.linspace(self.simulation_options["dispersion_parameters"][i]+"_lower", self.simulation_options["dispersion_parameters"][i]+"_upper", self.simulation_options["dispersion_bins"][i]))
                    weight_arrays.append([1/self.simulation_options["dispersion_bins"][i]]*self.simulation_options["dispersion_bins"][i])
            elif self.simulation_options["dispersion_distributions"][i]=="normal":
                if (self.simulation_options["dispersion_parameters"][i]+"_mean" not in self.dim_dict) or (self.simulation_options["dispersion_parameters"][i]+"_std" not in self.dim_dict):
                    raise ValueError("Uniform distribution requires "+self.simulation_options["dispersion_parameters"][i]+"_mean and " + self.simulation_options["dispersion_parameters"][i]+"_std")
                else:
                    param_mean=nd_dict[self.simulation_options["dispersion_parameters"][i]+"_mean"]
                    param_std=self.simulation_options["dispersion_parameters"][i]+"_std"
                    min_val=norm.ppf(1e-11, loc=param_mean, scale=param_scale)
                    max_val=norm.ppf(1-1e-11, loc=param_mean, scale=param_scale)
                    param_vals=np.linspace(min_val, max_val, self.simulation_options["dispersion_bins"][i])
                    param_weights=np.zeros(self.simulation_options["dispersion_bins"][i])
                    param_weights[0]=norm.cdf(param_vals[0],loc=param_mean, scale=self.param_std)
                    for j in range(1, self.simulation_options["dispersion_bins"][i]):
                        param_weights[j]=norm.cdf(param_vals[j],loc=param_mean, scale=self.param_std)-norm.cdf(param_vals[j-1],loc=param_mean, scale=param_std)
                    weight_values.append(param_vals)
                    weight_arrays.append(param_weights)
            elif self.simulation_options["dispersion_distributions"][i]=="lognormal":
                if (self.simulation_options["dispersion_parameters"][i]+"_shape" not in self.dim_dict) or (self.simulation_options["dispersion_parameters"][i]+"_loc" not in self.dim_dict) or (self.simulation_options["dispersion_parameters"][i]+"_scale" not in self.dim_dict):
                    raise ValueError("Uniform distribution requires "+self.simulation_options["dispersion_parameters"][i]+"_shape and " + self.simulation_options["dispersion_parameters"][i]+"_loc and "  + self.simulation_options["dispersion_parameters"][i]+"_scale")
                else:
                    param_loc=nd_dict[self.simulation_options["dispersion_parameters"][i]+"_loc"]
                    param_shape=nd_dict[self.simulation_options["dispersion_parameters"][i]+"_shape"]
                    param_scale=nd_dict[self.simulation_options["dispersion_parameters"][i]+"_scale"]
                    min_val=lognorm.ppf(1e-11, param_shape, loc=param_loc, scale=param_scale)
                    max_val=lognorm.ppf(1-1e-11, param_shape, loc=param_loc, scale=param_scale)
                    param_vals=np.linspace(min_val, max_val, self.simulation_options["dispersion_bins"][i])
                    param_weights=np.zeros(self.simulation_options["dispersion_bins"][i])
                    param_weights[0]=lognorm.cdf(param_vals[0],param_shape, loc=param_loc, scale=param_scale)
                    for j in range(1, self.simulation_options["dispersion_bins"][i]):
                        param_weights[j]=norm.cdf(param_vals[j],param_shape, loc=param_loc, scale=param_scale)-norm.cdf(param_vals[j-1],param_shape, loc=param_loc, scale=param_scale)
                    weight_values.append(param_vals)
                    weight_arrays.append(param_weights)
            else:
                raise KeyError(self.simulation_options["dispersion_distributions"][i]+" distribution not implemented")





    def kinetic_dispersion(self):
        #print self.nd_param.k0_shape, self.nd_param.k0_loc, self.nd_param.k0_scale
        k0_weights=np.zeros(self.simulation_options["dispersion_bins"])
        k_start=lognorm.ppf(1e-9, self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        k_end=lognorm.ppf(1-1e-5, self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        k0_vals=np.linspace(k_start,k_end, self.simulation_options["dispersion_bins"])
        k0_weights[0]=lognorm.cdf(k0_vals[0], self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        for k in range(1, len(k0_weights)):
            k0_weights[k]=lognorm.cdf(k0_vals[k], self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)-lognorm.cdf(k0_vals[k-1], self.nd_param.k0_shape, 0, scale=self.nd_param.k0_scale)
        #plt.title("k0")
        return k0_vals, k0_weights

    def therm_dispersion(self):
        self.e0_min=norm.ppf(1e-11, loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        self.e0_max=norm.ppf(1-(1e-11), loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
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
    def alpha_dispersion(self):
        self.e0_min=norm.ppf(1e-11, loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
        self.e0_max=norm.ppf(1-(1e-11), loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
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
    def numerical_plots(self, solver):
        self.debug_time=self.simulation_options["numerical_debugging"]
        time_series=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], self.debug_time, self.bounds_val)
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

    def simulate(self,parameters, frequencies):
        if len(parameters)!= len(self.optim_list):
            print(self.optim_list)
            print(parameters)
            raise ValueError('Wrong number of parameters')

        if self.simulation_options["label"]=="cmaes":
            normed_params=self.change_norm_group(parameters, "un_norm")
        else:
            normed_params=copy.deepcopy(parameters)
        #print(list(normed_params))
        for i in range(0, len(self.optim_list)):
            self.dim_dict[self.optim_list[i]]=normed_params[i]
        if "phase_only" in self.simulation_options and self.simulation_options["phase_only"]==True:
            self.dim_dict["cap_phase"]=self.dim_dict["phase"]
        self.nd_param=params(self.dim_dict)
        self.nd_param_dict=self.nd_param.nd_param_dict
        if self.simulation_options["numerical_method"]=="Brent minimisation":
            solver=isolver_martin_brent.brent_current_solver
        elif self.simulation_options["numerical_method"]=="Newton-Raphson":
            solver=isolver_martin_NR.NR_current_solver
            if self.simulation_options["method"]=="dcv":
                raise ValueError("Newton-Raphson dcv simulation not implemented")
        else:
            raise ValueError('Numerical method not defined')
        start=time.time()
        if self.simulation_options["numerical_debugging"]!=False:
            self.numerical_plots(solver)
        else:
            if self.simulation_options["dispersion"]==True:
                #print("dispersion")
                time_series=self.paralell_disperse(solver)
            else:
                time_series=solver(self.nd_param_dict, self.time_vec, self.simulation_options["method"],-1, self.bounds_val)
        #self.time_array.append(time.time()-start)
        if self.simulation_options["no_transient"]!=False:
            time_series=time_series[self.time_idx:]
        time_series=np.array(time_series)

        #print(time.time()-start, len(time_series))
        #time_series=self.define_voltages()
        if self.simulation_options["likelihood"]=='fourier':
            filtered=self.kaiser_filter(time_series)
            if (self.simulation_options["test"]==True):
                print(list(normed_params))
                plt.plot(self.secret_data_fourier, label="data")
                plt.plot(filtered , alpha=0.7, label="numerical")
                plt.show()

            return filtered
        elif self.simulation_options["likelihood"]=='timeseries':
            if self.simulation_options["test"]==True:
                print(list(normed_params))
                if self.simulation_options["experimental_fitting"]==True:
                    plt.subplot(1,2,1)
                    plt.plot(self.other_values["experiment_voltage"],time_series)
                    plt.plot(self.other_values["experiment_voltage"],self.secret_data_time_series, alpha=0.7)
                    plt.subplot(1,2,2)
                    plt.plot(self.other_values["experiment_time"],time_series)
                    plt.plot(self.other_values["experiment_time"],self.secret_data_time_series, alpha=0.7)
                    plt.show()
                else:
                    plt.plot(self.time_vec[self.time_idx:], time_series)
                    plt.plot(self.time_vec[self.time_idx:], self.secret_data_time_series)
                    plt.show()
            return (time_series)
class paralell_class:
    def __init__(self, params, times, method, bounds, solver):
        self.params=params
        self.times=times
        self.method=method
        self.bounds=bounds
        self.solver=solver
    def paralell_simulate(self, weight_val_entry):
        start=time.time()
        self.sim_params=copy.deepcopy(self.params)
        for i in range(len(weight_val_entry[0])):
            self.sim_params[weight_val_entry[0][i]]=weight_val_entry[1][i]
        time_series=self.solver(self.sim_params, self.times, self.method,-1, self.bounds)
        time_series=np.multiply(time_series, weight_val_entry[2])
        return (time_series)
    def paralell_dispersion(self, weight_list):
        p = mp.Pool(4)
        start1=time.time()
        sc = p.map_async(self,  [weight for weight in weight_list])
        start=time.time()
        results=sc.get()
        p.close()
        disped_time=np.sum(results, axis=0)
        start2=time.time()
        return disped_time
    def __call__(self, x):
        return self.paralell_simulate(x)
