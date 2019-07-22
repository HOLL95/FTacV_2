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
    def empty(self, arg):
        return arg
    def plot_harmonics(self, times, method, **kwargs):
        if method=="abs":
            a=abs
        else:
            a=self.empty
        harmonics_list=[]
        harmonics_labels=[]
        for key, value in kwargs.items():
            harmonics_list.append(value)
            harmonics_labels.append(key)
        fig, ax=plt.subplots(self.num_harmonics,1)
        for i in range(0, self.num_harmonics):
            for j in range(0, len(harmonics_labels)):
                    ax[i].plot(times, a(harmonics_list[j][i,:]), label=harmonics_labels[j])
        ax[i].yaxis.set_label_position("right")
        ax[i].set_ylabel(str(self.harmonics[i]), rotation=0)
        plt.legend()
        plt.show()
    def harmonics_and_time(self, times,title, method, **kwargs):
        names=kwargs.keys()
        titles=[]
        for i in range(0, len(names)):
            if ("_" in names[i]) and ("harmonics" in names[i]):
                titles.append(names[i][:names[i].index("_")])
        fig=plt.figure(num=None, figsize=(15, 6), dpi=80, facecolor='w', edgecolor='k')
        if method=="abs":
            a=abs
        else:
            a=self.empty
        time_list=[]
        time_labels=[]
        harmonics_list=[]
        harmonics_labels=[]
        for key, value in kwargs.items():
            if "harmonics" in key:
                harmonics_list.append(value)
                harmonics_labels.append(key)
            elif "alpha" in key:
                alpha=value
            else:
                time_list.append(value)
                time_labels.append(key)

        harm_axes=[]
        harm_len=2
        fig.text(0.03, 0.5, 'Current(A)', ha='center', va='center', rotation='vertical')

        for i in range(0,self.num_harmonics):
            harm_axes.append(plt.subplot2grid((self.num_harmonics,harm_len*2), (i,0), colspan=harm_len))
            for j in range(0, len(titles)):
                idx=harmonics_labels.index(titles[j]+"_harmonics")
                harm_axes[i].plot(times, a(harmonics_list[idx][i,:]), label=titles[j])
            harm_axes[i].yaxis.set_label_position("right")
            harm_axes[i].set_ylabel(str(self.harmonics[i]), rotation=0)
        harm_axes[i].legend()
        harm_axes[i].set_xlabel("Time(s)")

        time_ax=plt.subplot2grid((self.num_harmonics,harm_len*2), (0,harm_len), rowspan=self.num_harmonics, colspan=harm_len)
        for j in range(0, len(titles)):
            idx=time_labels.index(titles[j]+"_time_series")
            time_ax.plot(times, time_list[idx], label=titles[j], alpha=0.7)
        time_ax.set_ylabel("Current(A)")
        time_ax.set_xlabel("Time(s)")
        plt.legend()
        plt.suptitle(title)
        plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92, wspace=0.23)
        plt.show()
    def harmonics_and_voltages(self, times,voltages, title, method, **kwargs):
        names=kwargs.keys()
        titles=[]
        for i in range(0, len(names)):
            if ("_" in names[i]) and ("harmonics" in names[i]):
                titles.append(names[i][:names[i].index("_")])
        fig=plt.figure(num=None, figsize=(15, 6), dpi=80, facecolor='w', edgecolor='k')
        if method=="abs":
            a=abs
        else:
            a=self.empty
        voltage_list=[]
        voltage_labels=[]
        harmonics_list=[]
        harmonics_labels=[]
        for key, value in kwargs.items():
            if "harmonics" in key:
                harmonics_list.append(value)
                harmonics_labels.append(key)

            else:
                voltage_list.append(value)
                voltage_labels.append(key)
        print voltage_labels
        harm_axes=[]
        harm_len=2
        fig.text(0.03, 0.5, 'Current(A)', ha='center', va='center', rotation='vertical')
        for i in range(0,self.num_harmonics):
            harm_axes.append(plt.subplot2grid((self.num_harmonics,harm_len*2), (i,0), colspan=harm_len))
            for j in range(0, len(titles)):
                idx=harmonics_labels.index(titles[j]+"_harmonics")
                harm_axes[i].plot(times, a(harmonics_list[idx][i,:]), label=titles[j])
            harm_axes[i].yaxis.set_label_position("right")
            harm_axes[i].set_ylabel(str(self.harmonics[i]), rotation=0)
        harm_axes[i].legend()
        harm_axes[i].set_xlabel("Time(s)")
        time_ax=plt.subplot2grid((self.num_harmonics,harm_len*2), (0,harm_len), rowspan=self.num_harmonics, colspan=harm_len)
        for j in range(0, len(voltage_list)):
            print titles[j]
            idx=voltage_labels.index(titles[j]+"_time_series")
            if titles[j].lower()=="experimental":
                time_ax.plot(voltages, voltage_list[idx], label=titles[j], alpha=0.7)
            else:
                time_ax.plot(voltages, voltage_list[idx], label=titles[j])

        time_ax.set_ylabel("Current(A)")
        time_ax.set_xlabel("Potential(V)")
        plt.legend()
        plt.suptitle(title)
        plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92, wspace=0.23)
        plt.show()
