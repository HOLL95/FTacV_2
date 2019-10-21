import numpy as np
import os
import math
class FTACV_initialisation:
    def __init__(self, experimental_fitting, **kwargs):
        if experimental_fitting==True:
            self.voltage_results, self.current_results, self.time_results=self.file_opener(kwargs["file_dict"], kwargs["dec_amount"])

        self.generic_params()
    def file_opener(self, file_dict, dec_amount):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        data_path="/Experiment_data"
        folder="/Carbon"
        path=dir_path+data_path+folder
        files= os.listdir(path)
        voltage_results={}
        time_results={}
        current_results={}
        experiments=list(file_dict.keys())
        for data in files:
            for keys in experiments:
                Method=keys
                type=file_dict[keys]
                if (Method in data)  and (type[0] in data):
                    file=path+"/"+data
                    results=np.loadtxt(file)
                    current_results[Method]=results[0::dec_amount, 1]
                    time_results[Method]=results[0::dec_amount, 0]#
                    break
                elif (Method in data)  and (type[1] in data):
                    file=path+"/"+data
                    results=np.loadtxt(file)
                    voltage_results[Method]=results[0::dec_amount, 1]#
                    break
        return voltage_results, current_results, time_results
    def generic_params(self):
        de=300e-3
        estart=260e-3-de
        ereverse=estart+2*de
        self.generic_noramp_params={
            "E_0":0.25,
            'E_start': estart, #(starting dc voltage - V)
            'E_reverse': ereverse,
            'omega':5,#8.88480830076,  #    (frequency Hz)
            'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
            'v': 10.36e-3,   #       (scan rate s^-1)
            'area': 0.07, #(electrode surface area cm^2)
            'Ru': 1.0,  #     (uncompensated resistance ohms)
            'Cdl': 1e-5, #(capacitance parameters)
            'CdlE1': 0,#0.000653657774506,
            'CdlE2': 0,#0.000245772700637,
            'CdlE3': 0,#1.10053945995e-06,
            'gamma': 1e-10,
            "original_gamma":1e-10,        # (surface coverage per unit area)
            'k_0': 1.0, #(reaction rate s-1)
            'alpha': 0.5,
            "E0_mean":0.2,
            "E0_std": 0.09,
            "k0_shape":0.954,
            "k0_loc":100,
            "k0_scale":50,
            "k0_range":1e3,
            "cap_phase":3*(math.pi/2),
            'sampling_freq' : (1.0/200),
            'phase' : 3*(math.pi/2),
            "time_end": None,
            'num_peaks': 10
        }
        self.generic_ramped_params={
            "E_0":0.25,
            'E_start': -180e-3, #(starting dc voltage - V)
            'E_reverse': 620e-3,    #  (reverse dc voltage - V)
            'omega':5,#8.88480830076,  #    (frequency Hz)
            'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
            'v': 29.8e-3,   #       (scan rate s^-1)
            'area': 0.07, #(electrode surface area cm^2)
            'Ru': 1.0,  #     (uncompensated resistance ohms)
            'Cdl': 1e-5, #(capacitance parameters)
            'CdlE1': 0,#0.000653657774506,
            'CdlE2': 0,#0.000245772700637,
            'CdlE3': 0,#1.10053945995e-06,
            'gamma': 1e-10,
            "original_gamma":1e-10,          # (surface coverage per unit area)
            'k_0': 1.0, #(reaction rate s-1)
            "E0_mean":0.2,
            "E0_std": 0.09,
            "k0_shape":0.954,
            "k0_loc":100,
            "k0_scale":50,
            "k0_range":1e3,
            'alpha': 0.5,
            'sampling_freq' : (1.0/200),
            "cap_phase":0,
            'phase' : 0,
            'time_end':1000,
            'num_peaks': 600
        }
