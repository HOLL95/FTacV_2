#!/usr/bin/env python
import math
import copy
class params:
    def e0(self, value, flag):
        if flag=='re_dim':
            self.E_0=value*self.c_E0
        elif flag == 'non_dim':
            self.E_0=value/self.c_E0
    def k0(self, value, flag):
        if flag=='re_dim':
            self.k_0=value/self.c_T0
        elif flag == 'non_dim':
            self.k_0=value*self.c_T0
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
            #print "Changing the sampling freq to " + str(self.sampling_freq)

    def __init__(self,param_dict):
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
        self.T=(273+25)
        self.F=96485.3328959
        self.R=8.314459848
        self.c_E0=(self.R*self.T)/self.F
        self.c_T0=abs(self.c_E0/self.v)
        self.c_I0=(self.F*self.area*self.gamma)/self.c_T0
        self.c_Gamma=param_dict['gamma']
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
        keys=param_dict.keys()
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
