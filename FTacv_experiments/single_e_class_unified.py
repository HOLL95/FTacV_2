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
        if "GH_quadrature" in self.simulation_options:
            if self.simulation_options["GH_quadrature"]==True:
                self.GH_nodes, self.GH_weights=np.polynomial.hermite.hermgauss(self.simulation_options["dispersion_bins"])
                self.normal_GH_weights=[(1/math.sqrt(math.pi))*weight for weight in self.GH_weights]
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
        disp_flags=["mean", "scale", "upper"]
        disp_check=[[y in x for y in disp_flags] for x in self.optim_list]
        if True in [True in x for x in disp_check]:
            self.simulation_options["dispersion"]=True
            distribution_flags=["normal", "lognormal", "uniform"]
            self.simulation_options["dispersion_parameters"]=[]
            self.simulation_options["dispersion_distributions"]=[]
            for i in range(0, len(self.optim_list)):
                count=0
                for j in range(0, len(disp_flags)):
                    if count>1:
                        raise ValueError("Multiple dispersion flags in "+self.optim_list[i])
                    if disp_flags[j] in self.optim_list[i]:
                        index=self.optim_list[i].find("_"+disp_flags[j])
                        self.simulation_options["dispersion_parameters"].append(self.optim_list[i][:index])
                        self.simulation_options["dispersion_distributions"].append(distribution_flags[j])
                        count+=1
        else:
            self.simulation_options["dispersion"]=False
            #if type(self.simulation_options["dispersion_bins"]) is not list:
            #    self.simulation_options["dispersion_bins"]=[self.simulation_options["dispersion_bins"]]*len(self.simulation_options["dispersion_parameters"])

        if "phase" in optim_list and "cap_phase" not in optim_list:
            self.simulation_options["phase_only"]=True
        else:
            self.simulation_options["phase_only"]=False
        if "alpha_mean" in optim_list or "alpha_std" in optim_list:
            self.simulation_options["alpha_dispersion"]="normal"
        else:
            if "alpha_dispersion" in self.simulation_options:
                del self.simulation_options["alpha_dispersion"]
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
    def kaiser_filter(self, time_series, harmonical=False):
        frequencies=self.frequencies
        L=len(time_series)
        window=np.hanning(L)
        time_series=np.multiply(time_series, window)
        f=np.fft.fftfreq(len(time_series), self.time_vec[1]-self.time_vec[0])
        Y=np.fft.fft(time_series)
        #Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        top_hat=copy.deepcopy(Y[0:len(frequencies)])
        scale_flag=False
        true_harm=self.nd_param.omega*self.nd_param.c_T0

        if "fourier_scaling" in self.simulation_options:
            print(self.simulation_options["fourier_scaling"])
            if self.simulation_options["fourier_scaling"]!=None:
                scale_flag=True
        if sum(np.diff(self.harmonic_range))!=len(self.harmonic_range)-1 or scale_flag==True:
            results=np.zeros(len(top_hat), dtype=complex)
            for i in range(0, self.num_harmonics):
                true_harm_n=true_harm*self.harmonic_range[i]
                index=tuple(np.where((frequencies<(true_harm_n+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm_n-(self.nd_param.omega*self.filter_val))))
                if scale_flag==True:
                    filter_bit=abs(top_hat[index])
                    min_f=min(filter_bit)
                    max_f=max(filter_bit)
                    filter_bit=[self.normalise(x, [min_f, max_f]) for x in filter_bit]
                else:
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
        return abs(results)
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

    def normal_gh_transform(self,location, scale):
        try:
            if len(self.GH_nodes)!=self.simulation_options["dispersion_bins"]:
                self.GH_nodes, self.GH_weights=np.polynomial.hermite.hermgauss(self.simulation_options["dispersion_bins"])
                self.normal_GH_weights=[(1/math.sqrt(math.pi))*weight for weight in self.GH_weights]
            nodes=[(scale*math.sqrt(2)*node)+location for node in self.GH_nodes]
        except:
            self.GH_nodes, self.GH_weights=np.polynomial.hermite.hermgauss(self.simulation_options["dispersion_bins"])
            self.normal_GH_weights=[(1/math.sqrt(math.pi))*weight for weight in self.GH_weights]
            nodes=[(scale*math.sqrt(2)*node)+location for node in self.GH_nodes]
        return nodes

    def test_vals(self, parameters, likelihood, test=False):
        orig_likelihood=self.simulation_options["likelihood"]
        orig_label=self.simulation_options["label"]
        orig_test=self.simulation_options["test"]
        self.simulation_options["likelihood"]=likelihood
        self.simulation_options["label"]="MCMC"
        self.simulation_options["test"]=test
        results=self.simulate(parameters, self.frequencies)
        if sum(results)==0:
            print(self.simulation_options["dispersion"], self.simulation_options["alpha_dispersion"])
            raise ValueError("Not simulated, check options")
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
            counter=0
            if "alpha_dispersion" in self.simulation_options:
                if self.simulation_options["alpha_dispersion"]=="uniform":
                    alpha_bins=5
                    alpha_vals=np.linspace(self.param_bounds["alpha"][0],self.param_bounds["alpha"][1], alpha_bins)
                    alpha_weights=[1/alpha_bins]*alpha_bins
                elif self.simulation_options["alpha_dispersion"]=="normal":
                    if "GH_quadrature" in self.simulation_options:
                        if self.simulation_options["GH_quadrature"]==True:
                            alpha_vals=self.normal_gh_transform(location=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)
                            alpha_weights=self.normal_GH_weights
                            alpha_bins=len(alpha_vals)
                    else:
                        alpha_bins=16
                        alpha_min=max(0, norm.ppf(1e-4, loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std))
                        alpha_max=min(1,norm.ppf(1-(1e-4), loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std))
                        alpha_weights=np.zeros(alpha_bins)
                        alpha_vals=np.linspace(alpha_min, alpha_max, alpha_bins)
                        alpha_weights[0]=norm.cdf(alpha_vals[0], loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)
                        alpha_mids=np.zeros(alpha_bins)
                        alpha_mids[0]=alpha_vals[0]
                        for i in range(1, len(alpha_weights)):
                            alpha_weights[i]=norm.cdf(alpha_vals[i],loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)-norm.cdf(alpha_vals[i-1],loc=self.nd_param.alpha_mean, scale=self.nd_param.alpha_std)
                            alpha_mids[i]=(alpha_vals[i]+alpha_vals[i-1])/2
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
                    #plt.plot(voltages, np.multiply(time_series_current, e0_disp[i]))
                    #print(e0_disp[i])
        elif ("k0_shape" in self.optim_list):
            k0_vals, k0_disp=self.kinetic_dispersion()
            for i in range(0, self.simulation_options["dispersion_bins"]):
                self.nd_param_dict["k_0"]=k0_vals[i]
                time_series_current=solver(self.nd_param_dict, self.time_vec,self.simulation_options["method"], -1, self.bounds_val)
                time_series=np.add(time_series, np.multiply(time_series_current, k0_disp[i]))
        #plt.plot(k0_vals, k0_disp)
        #plt.show()


        return time_series
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
        if "GH_quadrature" in self.simulation_options:
            if self.simulation_options["GH_quadrature"]==True:
                e0_vals=self.normal_gh_transform(location=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
                e0_weights=self.normal_GH_weights

        else:
            self.e0_min=norm.ppf(1e-4, loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
            self.e0_max=norm.ppf(1-(1e-4), loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
            e0_weights=np.zeros(self.simulation_options["dispersion_bins"])
            e0_vals=np.linspace(self.e0_min, self.e0_max, self.simulation_options["dispersion_bins"])
            e0_midpoints=np.zeros(self.simulation_options["dispersion_bins"])
            e0_weights[0]=norm.cdf(e0_vals[0], loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
            e0_midpoints[0]=(e0_vals[0]+self.e0_min)/2
            for i in range(1, len(e0_weights)):
                e0_weights[i]=norm.cdf(e0_vals[i],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)-norm.cdf(e0_vals[i-1],loc=self.nd_param.E0_mean, scale=self.nd_param.E0_std)
                e0_midpoints[i]=(e0_vals[i]+e0_vals[i-1])/2
            e0_vals=e0_midpoints
            #midpoint=len(e0_weights)//2
            #e0_weights[midpoint:]=e0_weights[:midpoint]
            #range1=np.arange(self.simulation_options["dispersion_bins"], 0, -1)
            #e0_weights=np.divide(range1, sum(range1))
            #print(e0_weights)
            #plt.plot(self.e_nondim(e0_vals), e0_weights)
            #plt.show()
        #for i in range(0, len(e0_midpoints)):
        #    plt.axvline(self.e_nondim(e0_midpoints[i]), color="black", linestyle="--")
        #plt.show()

        #print(list(e0_weights))
        #print(list(self.e_nondim(e0_vals)))
        #plt.axvline(self.e_nondim(self.nd_param.E0_mean))
        #plt.plot(self.e_nondim(e0_vals), e0_weights)
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
                if type(self.simulation_options["dispersion_bins"]) is list:
                    if len(self.simulation_options["dispersion_bins"])!=1:
                        raise ValueError("Not currently implemented Henry")
                    else:
                        self.simulation_options["dispersion_bins"]=self.simulation_options["dispersion_bins"][0]
                time_series=self.paralell_disperse(solver)
                #print(np.sum(np.subtract(time_series, time_series2)))
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
                print("HEY")
                print(list(normed_params))
                plt.plot(self.secret_data_fourier, label="data")
                plt.plot(filtered , alpha=0.7, label="numerical")
                plt.legend()
                plt.show()

            return filtered
        elif self.simulation_options["likelihood"]=='timeseries':
            if self.simulation_options["test"]==True:
                print(list(normed_params))
                if self.simulation_options["experimental_fitting"]==True:
                    plt.subplot(1,2,1)
                    plt.plot(self.other_values["experiment_voltage"],time_series)
                    #plt.plot(self.other_values["experiment_voltage"],self.secret_data_time_series, alpha=0.7)
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
