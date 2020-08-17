
import isolver
import math
import numpy as np
import matplotlib.pyplot as plt
from anal_params_class import params
import copy
import time
import scipy.special as bessel
class analytical_electron:
    def __init__(self, dim_paramater_dictionary, Ein, gamma_0):
        key_list=list(dim_paramater_dictionary.keys())
        self.nd_param=params(dim_paramater_dictionary)
        for i in range(0, len(key_list)):
            self.nd_param.non_dimensionalise(key_list[i], dim_paramater_dictionary[key_list[i]])
        nd_dict=self.nd_param.__dict__
        self.nu=((self.nd_param.E_reverse+self.nd_param.E_start)/2)-self.nd_param.E_0
        I_0_1_alpha=bessel.iv(0, (1-self.nd_param.alpha)*self.nd_param.d_E)
        I_0_alpha=bessel.iv(0, self.nd_param.alpha*self.nd_param.d_E)
        self.gamma_inf=I_0_1_alpha/(I_0_1_alpha+(np.exp(-self.nu)*I_0_alpha))
        self.sigma=(np.exp((1-self.nd_param.alpha)*self.nu)*I_0_1_alpha)+np.exp(-self.nd_param.alpha*self.nu)*I_0_alpha
        self.gamma_0=gamma_0
    def h(self, t):
        h=np.exp((1-self.nd_param.alpha)*self.nu)
        h=h*np.exp((1-self.nd_param.alpha)*self.nd_param.d_E*np.sin(self.nd_param.omega*t+self.nd_param.phase))
        return h
    def g(self, t):
        g=self.h(t)+np.exp(-self.nd_param.alpha*self.nu)*np.exp((-self.nd_param.alpha)*self.nd_param.d_E*np.sin(self.nd_param.omega*t+self.nd_param.phase))
        return g
    def i(self,t):
        i=self.h(t)-(self.gamma_inf+(self.gamma_0-self.gamma_inf)*np.exp(-self.sigma*t))*self.g(t)
        return i
    def gamma_t(self, t):
        gamma_t=self.gamma_inf+(self.gamma_0-self.gamma_inf)*np.exp(-self.sigma*t)
