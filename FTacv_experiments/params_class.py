#!/usr/bin/env python
import math
class params:
    def e0(self, value, flag):
        if flag=='re_dim':
            self.E_0=value*self.c_E0
            if self.dispersion==True:
                self.E0_std=self.E0_std*self.c_E0
                self.E0_mean=self.E0_mean*self.c_E0
        elif flag == 'non_dim':
            self.E_0=value/self.c_E0
            if self.dispersion==True:
                self.E0_std=self.E0_std/self.c_E0
                self.E0_mean=self.E0_mean/self.c_E0

    def k0(self, value, flag):
        if flag=='re_dim':
            self.k_0=value/self.c_T0
            if self.dispersion==True:
                self.k0_loc=self.k0_loc/self.c_T0
                self.k0_scale=self.k0_scale/self.c_T0
        elif flag == 'non_dim':
            self.k_0=value*self.c_T0
            if self.dispersion==True:
                self.k0_loc=self.k0_loc*self.c_T0
                self.k0_scale=self.k0_scale*self.c_T0
    def cdl(self, value, flag):
        if flag=='re_dim':
            self.Cdl=value*self.c_I0*self.c_T0/(self.area*self.c_E0)
        elif flag == 'non_dim':
            self.Cdl=value/self.c_I0/self.c_T0*(self.area*self.c_E0)
    def estart(self, value, flag):
        if flag=='re_dim':
            self.E_start=value*self.c_E0
        elif flag == 'non_dim':
            self.E_start=value/self.c_E0
    def erev(self, value, flag):
        if flag=='re_dim':
            self.E_reverse=value*self.c_E0
        elif flag == 'non_dim':
            self.E_reverse=value/self.c_E0
    def omega_d(self, value, flag):
        if flag=='re_dim':
            self.omega=value/(2*math.pi*self.c_T0)
        elif flag == 'non_dim':
            self.nd_omega=value*(2*math.pi*self.c_T0)

    def de(self, value, flag):
        if flag=='re_dim':
            self.d_E=value*self.c_E0
        elif flag == 'non_dim':
            self.d_E=value/self.c_E0
    def ru(self, value, flag):
        if flag=='re_dim':
            self.Ru=value*self.c_E0/self.c_I0
        elif flag == 'non_dim':
            self.Ru=value/self.c_E0*self.c_I0
    def Gamma(self, value, flag):
        if flag=='re_dim':
            self.gamma=value/(1/self.c_Gamma)
        elif flag == 'non_dim':
            self.gamma=value*(1/self.c_Gamma)
    def sf(self, value, flag):
        if flag=='re_dim':
            self.sampling_freq=value/((2*math.pi)/self.nd_omega)
        elif flag == 'non_dim':
            self.sampling_freq=value*((2*math.pi)/self.nd_omega)
    def V(self, value, flag):
        if flag=='re_dim':
            self.v=value
        elif flag == 'non_dim':
            self.v=value
    def Phase(self, value, flag):
        if flag=='re_dim':
            self.phase=value
        elif flag == 'non_dim':
            self.phase=value
    def Alpha(self, value, flag):
        if flag=='re_dim':
            self.alpha=value
        elif flag == 'non_dim':
            self.alpha=value
    def ce1(self, value, flag):
        if flag=='re_dim':
            self.CdlE1=value
        elif flag == 'non_dim':
            self.CdlE1=value
    def ce2(self, value, flag):
        if flag=='re_dim':
            self.CdlE2=value
        elif flag == 'non_dim':
            self.CdlE2=value
    def ce3(self, value, flag):
        if flag=='re_dim':
            self.CdlE3=value
        elif flag == 'non_dim':
            self.CdlE3=value
    def __init__(self,param_dict):
        if "E0_std" in param_dict:
            self.dispersion=True
        else:
            self.dispersion=False
        self.E_0=param_dict['E_0']
        self.k_0=param_dict['k_0']
        self.Cdl=param_dict['Cdl']
        self.CdlE1=param_dict['CdlE1']
        self.CdlE2=param_dict['CdlE2']
        self.CdlE3=param_dict['CdlE3']
        self.d_E=param_dict['d_E']
        self.E_start=param_dict['E_start']
        self.E_reverse=param_dict['E_reverse']
        self.omega=param_dict['omega']
        self.alpha=param_dict['alpha']
        self.Ru=param_dict['Ru']
        self.gamma=param_dict['gamma']
        self.v=param_dict['v']
        self.sampling_freq=param_dict['sampling_freq']
        self.area=param_dict['area']
        self.phase=param_dict['phase']
        self.num_peaks=param_dict['num_peaks']
        self.time_end=param_dict['time_end']
        self.cap_phase=param_dict["cap_phase"]
        self.c_Gamma=param_dict['original_gamma']
        if self.dispersion==True:
            self.E0_std=param_dict["E0_std"]
            self.E0_mean=param_dict["E0_mean"]
            self.k0_shape=param_dict["k0_shape"]
            self.k0_loc=param_dict["k0_loc"]
            self.k0_scale=param_dict["k0_scale"]
            self.k0_range=param_dict["k0_range"]

        self.T=(273+25)
        self.F=96485.3328959
        self.R=8.314459848
        self.c_E0=(self.R*self.T)/self.F
        self.c_T0=abs(self.c_E0/self.v)
        self.c_I0=(self.F*self.area*self.c_Gamma)/self.c_T0

        self.method_switch={
                            'e_0':self.e0,
                            'k_0':self.k0,
                            'cdl':self.cdl,
                            'e_start' :self.estart,
                            'e_reverse': self.erev,
                            'omega':self.omega_d,
                            'd_e' :self.de,
                            'ru':self.ru,
                            'gamma':self.Gamma,
                            'sampling_freq':self.sf
                            }
        keys=sorted(param_dict.keys())
        for i in range(0, len(keys)):
            if keys[i].lower() in self.method_switch:
                self.non_dimensionalise(keys[i], param_dict[keys[i]])


    def non_dimensionalise(self, name,name_value):
        function = self.method_switch[name.lower()]
        function(name_value, 'non_dim')
    def re_dimensionalise(self,name, name_value):
        if name.lower() in self.method_switch:
                function = self.method_switch[name.lower()]
                function(name_value, 're_dim')
        else:
            raise ValueError(name + " not in param list!")
