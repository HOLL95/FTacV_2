import  matplotlib.pyplot as plt
import pints.plot
import numpy as np
import copy
class harmonics:
    def __init__(self, harmonics, input_frequency, filter_val):
        self.harmonics=harmonics
        self.num_harmonics=len(harmonics)
        self.input_frequency=input_frequency
        self.filter_val=filter_val
        print "initialised!"
    def generate_harmonics(self, times, data):
        L=len(data)
        window=np.hanning(L)
        time_series=np.multiply(data, window)
        f=np.fft.fftfreq(len(time_series), times[1]-times[0])
        Y=np.fft.fft(time_series)
        #plt.plot(self.test_frequencies, Y[:len(self.test_frequencies)])
        #plt.plot(f, Y)
        #plt.show()
        #Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        last_harm=(self.harmonics[-1]*self.input_frequency)
        frequencies=f[np.where((f>0) & (f<(last_harm+(0.5*self.input_frequency))))]
        top_hat=(copy.deepcopy(Y[0:len(frequencies)]))
        harmonics=np.zeros((self.num_harmonics, len(time_series)), dtype="complex")
        for i in range(0, self.num_harmonics):
            true_harm=self.harmonics[i]*self.input_frequency
            filter_bit=top_hat[np.where((frequencies<(true_harm+(self.input_frequency*self.filter_val))) & (frequencies>true_harm-(self.input_frequency*self.filter_val)))]
            harmonics[i,0:len(filter_bit)]=filter_bit
            harmonics[i,:]=((np.fft.ifft(harmonics[i,:])))
        return harmonics

    def plot_harmonics(self, times, harmonics, harmonics2=[False, False], method="abs", label1="", label2=""):
        fig, ax=plt.subplots(self.num_harmonics,1)
        for i in range(0, self.num_harmonics):
            if method=="abs":
                ax[i].plot(times, abs(harmonics[i,:]), label=label1)
            else:
                ax[i].plot(times, (harmonics[i,:]),label=label1)
            if label1!="":
                plt.legend()
            if harmonics2.any()!=False:
                if method=="abs":
                    ax[i].plot(times, abs(harmonics2[i,:]),label=label2)
                else:
                    ax[i].plot(times, (harmonics2[i,:]),label=label2)
                if label2!="":
                    plt.legend()
            ax[i].yaxis.set_label_position("right")
            ax[i].set_ylabel(str(self.harmonics[i]), rotation=0)
        plt.show()
    def comparison_harmonics_plot(self, times, harmonics, harmonics2, axis):
        for i in range(0, self.num_harmonics):
            axis.plot(times, abs(np.subtract(harmonics, harmonics2)))
