import os
import sys
import numpy as np
import matplotlib.pyplot as plt
dir_path = os.path.dirname(os.path.realpath(__file__))
slash_idx=[i for i in range(len(dir_path)) if dir_path[i]=="/"]
one_above=dir_path[:slash_idx[-1]]
sys.path.insert(1, one_above)
from single_e_class_unified import single_electron
Electrode="Yellow"
path=("/").join([dir_path , Electrode])
files=os.listdir(path)#
file_numbers=[]
counter=1
for names in files:
    plt.subplot(1,3,counter)
    file_numbers.append(names[names.index(".")-1])
    result=single_electron(path+"/"+names)
    #for i in range(0, len(self.save_dict[params])):
    plt.plot(result.other_values["experiment_time"], result.other_values["experiment_current"])
    for i in range(0, 3):
        simulations=(result.saved_param_simulate(result.save_dict["params"][i]))
        plt.plot(result.other_values["experiment_time"], simulations)
        plt.plot(result.other_values["experiment_time"], np.subtract(simulations,result.other_values["experiment_current"]) )
    counter+=1
plt.show()
