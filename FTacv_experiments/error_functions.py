import numpy as np
import matplotlib.pyplot as plt
class errors:
    def __init__(self, data):
        self.data=data
    def RMSE(self, prediction, absol=False):
        self.ninv=1.0/len(self.data)
        if absol==False:
            return np.sqrt(self.ninv * np.sum(
                (prediction - self.data)**2))
        else:
            return np.sqrt(self.ninv * np.sum(
                abs(prediction - self.data)**2))
    def SSE(self, prediction, absol=False):
        if absol==False:
            return np.sum((((prediction - self.data)**2)),
                              axis=0)
        else:
            return np.sum(((np.abs(prediction - self.data)**2)),
                              axis=0)
    def RE(self, prediction):
        prediction=prediction[np.where(prediction!=0)]
        data=self.data[np.where(self.data!=0)]
        return np.sum(((abs(1 - (prediction/data))**2)),
                          axis=0)
