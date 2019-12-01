import numpy as np
import matplotlib.pyplot as plt
import math
import time
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
import os
import pickle
import pints
dir_path = os.path.dirname(os.path.realpath(__file__))
types=["current", "voltage"]
exp="Experimental-120919"
blanks=["Blank-131119", "Blank-110919"]
exps=["Varying_salt", "Yellow"]
experiments=["Ramped", "Noramp"]
def binary_file_reader(filename):
    file=open(filename, "rb")
    binary_array=[]
    for line in file:
        string=line.decode("latin1")
        binary_array.append([float(string[0:string.index("\t")]), float(string[string.index("\t")+len("\t"):string.index("\r")])])
    return np.array(binary_array)
for j in range(0, len(experiments)):
    for i in range(0, 2):
        exp_type=blanks[i]
        #if exp_type==bla:
    #        extra="Blank/"
        #else:
        #    extra=""
        data_path="experiment_data_2/"+exp_type
        Experiment=exps[i]
        folder=experiments[j]


        type="current"
        type2="voltage"

        path=("/").join([dir_path, data_path, folder, Experiment])
        files= os.listdir(path)
        desired_conc="1"
        binary_files={}
        for filename in files:
            for type in types:
                if type in filename and desired_conc in filename:
                    binary_files[type]=path+"/"+filename



        experiment_data={}
        for type in types:
            experiment_data[type]=binary_file_reader(binary_files[type])
        dec_amount=1
        current_results1=experiment_data["current"][0::dec_amount, 1]
        time_results1=experiment_data["current"][0::dec_amount,0]
        voltage_results1=experiment_data["voltage"][0::dec_amount,1]
        plt.subplot(1,2, j+1)

        plt.ylabel("Current(nA)")
        if folder=="Ramped":
            plt.plot(time_results1[:-50], np.multiply(current_results1[:-50], 1e9), alpha=0.5, label=blanks[i])
            plt.xlabel("Time(s)")
        else:
            plt.plot(voltage_results1[:-500], np.multiply(current_results1[:-500], 1e9),label=blanks[i])
            plt.xlabel("Voltage(V)")
        plt.legend()
plt.show()
