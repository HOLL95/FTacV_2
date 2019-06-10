import math
import numpy as np
import matplotlib.pyplot as plt
from params_class import params
from single_e_class_noramp import single_electron
import copy
import time
class electrochem_funcs:
    def __init__(self, nd_param_class ):
        self.nd_param=nd_param_class
    def kox(self, potential, current):
        Er=potential-(self.nd_param.Ru*current)
        exp11=np.exp(-self.nd_param.alpha*(Er-self.nd_param.E_0))
        return self.nd_param.k_0*exp11
    def kred(self, potential, current):
        Er=potential-(self.nd_param.Ru*current)
        exp12=np.exp((1-self.nd_param.alpha)*(Er-self.nd_param.E_0))
        return self.nd_param.k_0*exp12
