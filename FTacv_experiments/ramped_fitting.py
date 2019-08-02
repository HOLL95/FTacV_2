import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_ramped
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from single_e_class_ramped  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]

dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Carbon"
Method ="GC4_1_ramp"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        voltages=np.loadtxt(path+"/"+data)
dec_amount=64
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=voltages[0::dec_amount, 1]
Method ="GC4_1_cv"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        voltages=np.loadtxt(path+"/"+data)
dec_amount=64
current_results2=results[0::dec_amount, 1]
time_results2=results[0::dec_amount, 0]
voltage_results2=voltages[0::dec_amount, 1]
plt.plot(voltage_results, current_results)
plt.plot(voltage_results2, current_results2)
plt.axvline(0.21)
plt.axvline(0.267)
plt.show()



param_list={
    "E_0":0.2,
    'E_start': -180e-3, #(starting dc voltage - V)
    'E_reverse': 620e-3,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 29.8e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 20.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-6, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,
    "original_gamma":1e-10,          # (surface coverage per unit area)
    'k_0': 100.0, #(reaction rate s-1)
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
    'num_peaks': 50
}
likelihood_options=["timeseries", "fourier"]
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "test": False,
    "likelihood":likelihood_options[1],
    "dispersion":False,
    "dispersion_bins":16,
    "label": "cmaes",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(2,7,1),
    "experiment_time": time_results,
    "experiment_current": current_results,
    "experiment_voltage":voltage_results,
    "bounds_val":2000,
    "signal_length":len(time_results),
}


ramp_fit=single_electron(param_list, simulation_options, other_values)
time_results=ramp_fit.other_values["experiment_time"]
current_results=ramp_fit.other_values["experiment_current"]
voltage_results=ramp_fit.other_values["experiment_voltage"]
likelihood_func=ramp_fit.kaiser_filter(current_results)
ramp_fit.pass_extra_data(current_results, likelihood_func)
ramp_fit.optim_list=["gamma"]
ramp_fit.simulation_options["dispersion"]=False
#normal=ramp_fit.test_vals([1e-10], likelihood="timeseries", test=False)
#ramp_fit.simulation_options["dispersion"]=True
#disped=ramp_fit.test_vals([1e-10], likelihood="timeseries", test=False)
#plt.plot(normal)
#plt.plot(disped)
#plt.show()
ramp_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl','gamma', 'omega', 'phase', 'alpha']
param_bounds={
    'E_0':[0.15, 0.3],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 6e2],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-2, 2],#0.000653657774506,
    'CdlE2': [-0.2,0.2],#0.000245772700637,
    'CdlE3': [0,0.1],#1.10053945995e-06,
    'gamma': [1e-11,1e-9],
    'k_0': [0, 5e2], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    'k_0': [0, 1e3], #(reaction rate s-1)
    "E0_mean":[0.0, 0.5],
    "E0_std": [0.01, 0.15],
    "k0_shape":[0,5],
    "k0_loc":[1, 1e4],
    "k0_scale":[0,2e3],
    "k0_range":[1e2, 1e4],
    'phase' : [0, 2*math.pi]
}
ramp_fit.optim_list=[]

harm_class=harmonics(other_values["harmonic_range"], ramp_fit.nd_param.omega*ramp_fit.nd_param.c_T0, 0.03)
ramp_fit.optim_list=["E0_mean", "E0_std","k_0", 'Ru', 'Cdl',"CdlE1","CdlE2",'gamma', 'omega', 'phase', 'alpha']


ramp_means_black=[0.20504878429957712, 0.07357984171424121, 126.17588866800264, 24.24736830211999, 9.999999999876686e-07, 0.11886527586092166, 0,1.3175840584818795e-10, 8.959496500290681, 6.2831853071795845, 0.599999999999731]
ramp_means_black_2=[0.21200070197828386, 0.06888959981526988, 133.96563492653507, 40.08177557223102, 3.226207450320691e-06, -0.021487125154184827, 0.0017931237883237632, 1.2669148148700962e-10, 8.959483328711777, 6.283185307173828, 0.7999999999999661]
ramp_means_carbon=[0.23192913053278295, 0.07597303082211063,133.999986524228, 20, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 2.8654939556021e-10, 8.959351751379364, 6.077678909557666, 0.7999999909807196]
ramp_free_means_black=[0.20504878429957712, 0.04692985835905884, 773.0039074468887, 1.1172494386860095e-06, 6.253912022948444e-06, 0.46590284463560927, -0.020672008906663236, 9.665438282298918e-11, 8.94055300929529,6.2831853071795845,  0.10000008017506497]
ramp_free_means_black_2=[0.22026089333976873, 0.04776183826475387, 1226.0003897156193, 2.6091841820962103e-10, 6.311657164346574e-06, 0.4861839637949368, -0.021712454485320096, 9.470125832724226e-11, 8.959351751379364, 6.2831853071795845, 0.10000000178756306]
ramp_free_means_carbon=[0.2445537156141517, 0.07597303082211063,100.9519493198850647, 1.6017819143290766e-07, 3.315492244571699e-05, 0.08995451289980627, -0.003381307141896925, 1.625196327083214e-10, 8.959351751379364,  6.2831853071795845, 0.5415684133427566]
ramp_free_means_carbon=[0.24437320720249422,0.9641045582560264, 4.4556439040686615e-08, 3.3222199217955305e-05, 0.09074940228372252, -0.0034091889329414538, 1.6095034096343578e-10, 8.940780387771508, 6.2831853071795845, 6.2831853071795845, 0.5417098085713957]
ramp_free_means_carbon_1=[0.231, 0.07,133.99,22.68,3e-5,-0.191,0.0128,2.865e-10,8.959351751379364, 0, 0.5]
ramp_free_means_carbon_1=[0.2255608974358974, 0.05391826923076923, 48.076923076923094, 163.46153846153842, 3.661858974358974e-05, -0.0064102564102563875, 0.0009615384615384581, 1.512019230769231e-10, 8.959351751379364, 0.0, 0.5]
ramp_free_means_carbon_1=[0.21834935897435898, 0.04971153846153847, 97.75641025641022, 598.5576923076923, 3.661858974358974e-05, -0.0064102564102563875, 0.0009615384615384581, 2.0673076923076927e-10, 8.959351751379364, 0.010069207223044163, 0.5]
ramp_free_means_carbon_1=[0.21834935897435898, 0.04971153846153847, 50.75641025641022, 50.5576923076923, 3.661858974358974e-05, -0.0064102564102563875, 0.0009615384615384581, 1.512019230769231e-10, 8.959351751379364, 0.010069207223044163, 0.5]
ramp_free_means_carbon_1=[0.23264569588050996, 0.0674603084243984, 84.4445611107851, 227.55356998656478, 1.9282854937765275e-08, 1.9129110746656064, 0.16331952396542287, 1.9030623303683534e-10, 8.959299190871032, 5.418246137182846, 0.40000000005905]


#ramp_free_means_carbon_1=[0.2824519230769231, 0.08087339743589742, 39.2628205128205, 140.3846153846154, 3.8060897435897424e-05, -0.191, 0.004487179487179471, 3.503125000000001e-10, 8.959351751379364, 0.3473876491950252, 0.5]

#ramp_free_means_carbon_1=[0.19911858974358976, 0.05949519230769229, 36.05769230769232, 44.71153846153848, 8.253205128205123e-06, 0.18910256410256432, -0.00801282051282054, 1.0281250000000005e-10, 8.959351751379364, 0.0, 0.5]
#THIS IS BLACK
#[0.21915064102564102, 0.05949519230769229, 36.05769230769232, 44.71153846153848, 8.253205128205123e-06, 0.18910256410256432, -0.00801282051282054, 1.0281250000000005e-10, 8.959351751379364, 0.0, 0.5]
#CARBON 50 dispersion_bins
#[0.234375, 0.06483974358974359, 63.301282051282044, 40.38461538461539, 3.8060897435897424e-05, -0.191, 0.004487179487179471, 1.512019230769231e-10, 8.959351751379364, 0.0, 0.5]

#ramp_free_means_carbon_1=[0.22916666666666669, 0.05861513036242604, 52.083333333333314, 0.65606508875739, 3e-05, -0.14423076923077005, 0.004006410256410242, 1.408894230769231e-10, 8.959351751379364, 0.0, 0.5]
#ramp_free_means_carbon_1=[0.231169872, 0.0697195513, 127.403846, 410.096154, 2.7884615384615376e-05, -0.00801282051282115, 0.001, 3.44759615e-10, 8.959483328711777, 0.0, 0.516666667]

#ramp_free_means_carbon_1=[0.26689122715624614,0.07, 97.72540069109925, 601.2211952415878, 3.379513907637756e-05, 0.10070581762902808, -0.0038788288927744663, 1.4557809348768636e-10, 8.940779258804525, 0,0.5]
#ramp_free_means_carbon_1=[0.26689122715624614,0.07, 97.72540069109925, 100.2211952415878, 3.379513907637756e-05, 0.10070581762902808*0, -0.0038788288927744663*0, 1.4557809348768636e-10, 8.940779258804525, 0, 0.5]
ramp_free_means_carbon_2=[0.2712720627064147, 53.00972798907347*100, 126.43581153952566, 3.160962997490096e-05*0, 0.14647947194537103*0, -0.005848319334033306*0, 1.0072445202882476e-10, 8.940709552792356,0, 0.5488793068724522]
ramp_free_means_carbon_2=[0.2528081478002836, 217.56996476308117, 150.7297882492157, 3.0166883670561975e-05, 0.14790906812787918, -0.005389008883307372, 1.1330727195554263e-10, 8.940664568706621,0, 0.40000000111206985]
cmaes_ramped_time=np.zeros(len(current_results))#ramp_fit.test_vals(ramp_means_carbon, likelihood="timeseries", test=False)
ramp_fit.simulation_options["dispersion"]=True
ramp_fit.optim_list=["E0_mean", "E0_std",'k_0', 'Ru',"Cdl","CdlE1", "CdlE2",'gamma', 'omega', "phase",  "alpha"]
#cmaes_rampfree_time=ramp_fit.test_vals(ramp_free_means_carbon_2, likelihood="timeseries", test=False)
#cmaes_rampfree_time_2=ramp_fit.test_vals(ramp_free_means_black_2, likelihood="timeseries", test=False)

data_harmonics=harm_class.generate_harmonics(time_results, current_results)

#ramp_harmonics=harm_class.generate_harmonics(time_results, cmaes_ramped_time)
#ramp_free_harmonics=harm_class.generate_harmonics(time_results, cmaes_rampfree_time)
#ramp_free_harmonics_2=harm_class.generate_harmonics(time_results, cmaes_rampfree_time_2)
#harm_class.plot_harmonics(time_results, method="abs", Experimental=data_harmonics, Ramped=ramp_harmonics)#, Ramp_free=ramp_free_harmonics)
#harm_class.harmonics_and_time(ramp_fit.t_nondim(time_results), folder, "abs", \
#                            Experimental_harmonics=ramp_fit.i_nondim(data_harmonics),Sinusoidal_harmonics=ramp_fit.i_nondim(ramp_free_harmonics),Ramped_harmonics=ramp_fit.i_nondim(ramp_harmonics),\
#                             Experimental_time_series=ramp_fit.i_nondim(current_results),Sinusoidal_time_series=ramp_fit.i_nondim(cmaes_rampfree_time), Ramped_time_series=ramp_fit.i_nondim(cmaes_ramped_time)) ##
test=ramp_fit.test_vals(ramp_free_means_carbon_1, likelihood="timeseries", test=False)

f_test=ramp_fit.kaiser_filter(test)
likelihood_func=ramp_fit.kaiser_filter(current_results)
print len(f_test), len(likelihood_func)
plt.plot(likelihood_func)
plt.plot(f_test)
plt.show()
error=np.power(np.subtract(f_test, likelihood_func),2)
error=np.sum(error)
print error
plt.plot(abs(np.fft.fft(current_results)[:len(current_results)/2+1]))
plt.plot(abs(np.fft.fft(test)[:len(test)/2+1]))
plt.show()
test_harmonics=harm_class.generate_harmonics(time_results, test)
harm_class.harmonics_and_time(ramp_fit.t_nondim(time_results), folder, "abs",\
                                Experimental_harmonics=ramp_fit.i_nondim(data_harmonics) , Experimental_time_series=ramp_fit.i_nondim(current_results),\
                                test_harmonics=ramp_fit.i_nondim(test_harmonics), test_time_series=ramp_fit.i_nondim(test))
#fig,ax=plt.subplots(1,1)

#l1,=(ax.plot(time_results, (test), lw=2))
#(ax.plot(time_results, (current_results)))v
fig, ax=plt.subplots(harm_class.num_harmonics+1, 1)
line_elements=[]
l1,=(ax[0].plot(voltage_results, test, lw=2))
ax[0].plot(voltage_results, current_results, lw=2, alpha=0.5)
line_elements.append(l1)
j=0
for i in range(1, harm_class.num_harmonics+1):
    l1,=(ax[i].plot(time_results, abs(test_harmonics[j,:]), lw=2))
    (ax[i].plot(time_results, abs(data_harmonics[j,:])))
    j+=1
    line_elements.append(l1)
true_data=ramp_free_means_carbon_1
"""
for i in range(0, len(ramp_fit.optim_list)):
    if ("Cdl" in ramp_fit.optim_list[i]) or ("gamma" in ramp_fit.optim_list[i]):
        continue
    else:
        param_bounds[ramp_fit.optim_list[i]][0]=true_data[i]*0.7
        param_bounds[ramp_fit.optim_list[i]][1]=true_data[i]*1.3
"""
axcolor = 'lightgoldenrodyellow'
slider_ax=[]
slider_ax_element=[]
for i in range(0, len(ramp_fit.optim_list)):
    slider_ax.append(plt.axes([0.25, 0.0+(i*0.02), 0.65, 0.01], facecolor=axcolor))
    slider_ax_element.append(Slider(slider_ax[i], ramp_fit.optim_list[i], param_bounds[ramp_fit.optim_list[i]][0], param_bounds[ramp_fit.optim_list[i]][1],ramp_free_means_carbon_1[i]) )



def update(val):
    params=np.zeros(len(ramp_fit.optim_list))
    for i in range(0, len(ramp_fit.optim_list)):
        params[i]=slider_ax_element[i].val
    print list(params)
    test=ramp_fit.test_vals(params, likelihood="timeseries", test=False)
    test_harmonics=harm_class.generate_harmonics(time_results, test)
    line_elements[0].set_ydata(test)
    j=0
    for i in range(1, harm_class.num_harmonics+1):
        line_elements[i].set_ydata(abs(test_harmonics[j,:]))
        j+=1
    fig.canvas.draw_idle()

for i in range(0, len(ramp_fit.optim_list)):
    slider_ax_element[i].on_changed(update)


resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')







plt.show()
ramp_fit.bounds_val=1e5


plt.plot(test*-1)
plt.plot(current_results)
plt.show()


ramp_fit.simulation_options["test"]=False
param_boundaries=np.zeros((2, ramp_fit.n_parameters()))
for i in range(0, ramp_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[ramp_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[ramp_fit.optim_list[i]][1]

ramp_fit.define_boundaries(param_boundaries)
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics,"abs", "numerical", "data")
#harm_class.harmonics_and_time(time_results, test_harmonics, test, data_time=current_results, harmonics2=data_harmonics,label1="numerical", label2="data", title="Black")
likelihood_func=ramp_fit.kaiser_filter(current_results)
dummy_times=np.linspace(0, 1, len(likelihood_func))
plt.plot(likelihood_func)
plt.show()
#cmaes_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
cmaes_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, likelihood_func)
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
x0=ramp_fit.change_norm_group(ramp_free_means_carbon_1, "norm")#abs(np.random.rand(ramp_fit.n_parameters()))

for i in range(0, 1):
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
#    x0=found_parameters
#    print found_parameters
print folder
cmaes_results=ramp_fit.change_norm_group(found_parameters, "un_norm")
print list(cmaes_results)
cmaes_time=ramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
print "TIME"
plt.plot(cmaes_time)
plt.plot(current_results)
plt.show()
cmaes_time=ramp_fit.test_vals(cmaes_results, likelihood="fourier", test=True)

"""
cmaes_time_prediction=ramp_fit.simulate(found_parameters,frequencies, "optimise", "timeseries", "yes" )
test_harmonics=harm_class.generate_harmonics(time_results, cmaes_time_prediction)
harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics,"abs", "numerical", "data")
#harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)
error=np.std(np.subtract(cmaes_prediction, likelihood_func))
#error=np.std(np.subtract(cmaes_time_prediction, current_results))
mcmc_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, likelihood_func)
#mcmc_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
updated_lb=np.append(cmaes_results*0.5, [0])#found_parameters[3]*0.97,
updated_ub=np.append(cmaes_results*1.5, [max(likelihood_func)])#found_parameters[3]*1.03,
updated_boundaries=[updated_lb, updated_ub]
updated_boundaries=np.sort(updated_boundaries, 0)
ramp_fit.define_boundaries(updated_boundaries)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                updated_boundaries[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(cmaes_results, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
ramp_fit.label="MCMC"
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()
pints.plot.trace(chains)
plt.show()
filename="ramped_results_2"
f=open(filename, "w")
np.save(f, chains)
f.close()
"""
