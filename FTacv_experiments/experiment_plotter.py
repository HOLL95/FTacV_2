import numpy as np
import os
import matplotlib.pyplot as plt
from harmonics_plotter import harmonics
import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import matplotlib.ticker as ticker
plt.rcParams.update({'font.size': 10})
harm_class=harmonics(list(range(1, 6)), 8.94, 0.05)
harmonic_no=4
max_harm=8
time_list=[]
omega=8.94
letters=['A', 'B', 'C', 'D', 'E', 'F', "G", 'H', 'I', 'J', 'K', 'L']
folder="/Experiment_data/Carbon"
names=["PostLinearRamp","GC4_1_ramp_cv", "GC4_1_cv"]
y_names=["Linear Sweep", "Ramped", "Sinusoidal"]
types=["Voltage input", "Current response"]
current="current"
voltage="voltage"
dcv_method="asA_4"
titles=["Voltage input", "Current response", "Absolute Fourier transform", "Harmonic "+str(harmonic_no) ]
def voltage_input(time, voltage, ax):
    ax.plot(time, voltage)
    ax.set_ylabel("Voltage(V)")
    ax.set_xlabel("Time(s)")
def time_series(voltage, data, ax):
    ax.plot(voltage, data*1000)
    ax.set_xlabel("Voltage(V)")
    ax.set_ylabel("Current(mA)")
def fourier(frequencies,freq_idx,  data, ax):
    hanning=np.ones(len(data))#np.hanning(len(data))
    f_data=np.fft.fft(np.multiply(data,hanning))
    log_f_data=(abs(f_data[freq_idx]))
    ax.plot(frequencies, log_f_data)
    ax.set_xlabel("Frequency(Hz)")
    ax.set_ylabel("Amplitude")
def harmonic_plot(time, harmonics, harmonic_no, ax):
    ax.plot(time, (harmonics*1000))
    #ax.title("Harmonic"+ str(harmonic_no))
    ax.set_xlabel("Time(s)")
    ax.set_ylabel("Current(mA)")

#ax.text(-0.1, 1.15, letters[letter_counter], transform=ax.transAxes,
    #fontsize=16, fontweight='bold', va='top', ha='right')
plt.rcParams.update({'font.size': 12})
fig, ax=plt.subplots(2, 2)
directory=os.path.dirname(os.path.realpath(__file__))
dec_amount=8
v_array=[]
c_array=[]
t_array=[]
for i in range(0, len(names)):
    if i ==0:
        path=directory+folder+"/dcv"
    else:
        path=directory+folder
    files =os.listdir(path)
    for filename in files:
        if (names[i] in filename)  and (current in filename):
            results=np.loadtxt(path+"/"+filename)
            current_results=results[::dec_amount, 1]
            time_results=results[::dec_amount, 0]
        elif (names[i] in filename)  and (voltage in filename):
            results_v=np.loadtxt(path+"/"+filename)
            voltage_results=results_v[::dec_amount, 1]
        elif (dcv_method in filename) and (names[i] in filename):
            dcv_results=open(path+"/"+filename)
            dcv_results=np.loadtxt(dcv_results, skiprows=1)
            time_results=dcv_results[:,0]
            voltage_results=dcv_results[:,1]
            current_results=dcv_results[:,2]

        if (i+1)==len(names):
            current_results=current_results[:10000]
            time_results=time_results[:10000]
            voltage_results=voltage_results[:10000]
    v_array.append(voltage_results)
    c_array.append(current_results)
    t_array.append(time_results)
    print(len(current_results), len(voltage_results), len(time_results))






for j in range(1, 3):
    for i in range(0, 2):
        voltage_results=v_array[j]
        current_results=c_array[j]
        time_results=t_array[j]
        axes=ax[j-1,i]
        if j==0:

            axes.set_title(types[i])
        if i==0:
            if j ==2:
                max=(1/8.94)*20
                idx=tuple(np.where(time_results<max))
                voltage_input(time_results[idx], voltage_results[idx], axes)
            else:
                voltage_input(time_results, voltage_results, axes)
            reduction=2465//3
            if (j+1)==3:
                axes.plot(time_results[:reduction], voltage_results[:reduction])
            pad = 5
            axes.annotate(y_names[j], xy=(0, 0.5), xytext=(-axes.yaxis.labelpad - pad, 0),
                xycoords=axes.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center', rotation="horizontal")

        elif i==1 and (j+1)==3:
            time_series(voltage_results[idx], current_results[idx], axes)
            axes.plot(voltage_results[:reduction], current_results[:reduction]*1000)
        else:
            time_series(voltage_results, current_results, axes)

        #elif j==2:
        #    frequencies=np.fft.fftfreq(len(current_results), time_results[1]-time_results[0])
        #    freq_idx=np.where((frequencies>0) & (frequencies<omega*12+(omega*0.5)))
        #    freqs=frequencies[freq_idx]
        #    fourier(freqs, freq_idx, current_results, axes)
        #elif j==3:
        #    harmonics=harm_class.generate_harmonics(time_results, current_results)
        #    harmonic_plot(time_results, abs(harmonics[harmonic_no, :]), harmonic_no, axes)






plt.show()
#plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92)
#fig.text(0.02, 0.21, 'Current(nA)', ha='center', va='center', rotation='vertical')
#fig.text(0.5, 0.21, 'Current(nA)', ha='center', va='center', rotation='vertical')
#plt.show()
