import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_noramp  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]
folder_options=["/Black", "/Red", "/Plain", "/Carbon", "/Gold/Large", "/Gold/Small"]
dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder=folder_options[3]#"/Black"#"/Gold/Large"#
Method ="GC4_1_cv"#"GoldLarge_1"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
print files
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        results2=np.loadtxt(path+"/"+data)
dec_amount=8
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=results2[0::dec_amount, 1]
#plt.plot(time_results, current_results)
#plt.show()

current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=results2[0::dec_amount, 1]
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-5*0, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,
    "original_gamma":1e-10,        # (surface coverage per unit area)
    'k_0': 10, #(reaction rate s-1)
    'alpha': 0.5,
    "E0_mean":0.2,
    "E0_std": 0.09,
    "k0_shape":0.954,
    "k0_loc":100,
    "k0_scale":50,
    "k0_range":1e3,
    "cap_phase":0,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
    "time_end": None,
    'num_peaks': 50
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
time_start=5/(param_list["omega"])
start_idx=np.where(time_results>time_start)
simulation_options={
    "no_transient":start_idx[0][0],
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":50,
    "test": False,
    "likelihood":likelihood_options[0],
    "numerical_method": solver_list[1],
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,9,1),
    "experiment_time": time_results,
    "experiment_current": current_results,
    "experiment_voltage":voltage_results,
    "bounds_val":20,
    "signal_length":int(2e4),
}
param_list['E_0']=2.36504746e-01#(param_list['E_reverse']-param_list['E_start'])/2
noramp_fit=single_electron(param_list, simulation_options, other_values)
time_results=noramp_fit.other_values["experiment_time"]
current_results=noramp_fit.other_values["experiment_current"]
voltage_results=noramp_fit.other_values["experiment_voltage"]
likelihood_func=noramp_fit.kaiser_filter(current_results)
test_voltages=noramp_fit.define_voltages()
noramp_fit.optim_list=[]
frequencies=noramp_fit.frequencies
noramp_fit.pass_extra_data(current_results, likelihood_func)
param_bounds={
    'E_0':[0.1, 0.4],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 1e3],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-2,2],#0.000653657774506,
    'CdlE2': [-0.1,0.1],#0.000245772700637,
    'CdlE3': [-0.05,0.05],#1.10053945995e-06,
    'gamma': [1e-11,4e-10],
    'k_0': [20, 1e4], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    "cap_phase":[0, 2*math.pi],
    "E0_mean":[0.0, 0.5],
    "E0_std": [0.01, 0.3],
    "k0_shape":[0,5],
    "k0_loc":[1, 1e4],
    "k0_scale":[0,2e3],
    "k0_range":[1e2, 1e4],
    'phase' : [0, 2*math.pi]
}
#'E_0', 'k_0' 'Ru', 'Cdl','gamma'
harm_class=harmonics(other_values["harmonic_range"], noramp_fit.nd_param.omega*noramp_fit.nd_param.c_T0, 0.1)
data_harmonics=harm_class.generate_harmonics(time_results, current_results)



noramp_fit.optim_list=['E0_mean', "E0_std",'k_0', 'Ru',"Cdl", "CdlE1", "CdlE2",'gamma', 'omega', "alpha", "phase", "cap_phase"]
#carbon_means=[0.24432865369586831, 0.03718591974278576, 5266.813117359024, 2.3351875006109577e-06, 3.904567925683219e-05, 0.020282592797868926,0, 2.886466778833776e-10, 8.940621635999642, 0.10000000291906254,]
carbon_means=[0.2629196040175071, 0.010000000426528552, 242.94186943883, 793.6555564492475, 4.020715996280745e-05, 0.012041277069243517,0, 1.702054870510513e-10, 8.94103486982429, 0.8999999976170034, 5.132595831851037, 4.446662466382718]
black_means=[0.20504878429957712, 0.04692985835905884, 773.0039074468887, 1.1172494386860095e-06, 6.253912022948444e-06, 0.46590284463560927, -0.020672008906663236, 9.665438282298918e-11, 8.94055300929529, 0.10000008017506497, 4.749140535109004, 4.029008370204711]
black_means_2=[0.22026089333976873, 0.04776183826475387, 1226.0003897156193, 2.6091841820962103e-10, 6.311657164346574e-06, 0.4861839637949368, -0.021712454485320096, 9.470125832724226e-11, 8.940448789535177, 0.10000000178756306, 4.742083569735882, 4.028465609498452]
black_means_3=[0.19813523414292655, 1.4877187307507795, 4.4556439040686615e-08, 1.1663407409020574e-06, 1.775725496685462, -0.07670261113515547, 1.2372296727209742e-10, 8.940758763885253, 5.622404458746398, 4.172942575228049, 0.5999999999990152]
carbon_means_2=[0.24437320720249422, 0.9641045582560264, 20, 3.3222199217955305e-05, 0.09074940228372252, -0.0034091889329414538, 1.6095034096343578e-10, 8.940780387771508, 5.698538006754413, 4.317679292947246, 0.5417098085713957]
carbon_means_3=[0.26689122715624614, 97.72540069109925, 601.2211952415878, 3.379513907637756e-05, 0.10070581762902808, -0.0038788288927744663, 1.4557809348768636e-10, 8.940779258804525, 5.149017610779459, 4.450212193055553]
ramp_means=[0.21200070197828386, 0.06888959981526988, 133.96563492653507, 40.08177557223102, 3.226207450320691e-06, -0.021487125154184827, 0.0017931237883237632, 1.2669148148700962e-10, 8.940448789535177, 0.7999999999999661,4.749140535109004, 4.028465609498452]
ramp_means_carbon_1=[0.24192913053278295, 9999.999986524228, 1.0428644292961376e-08, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 3.5654939556021e-10, 8.940621635999642, 4.749140535109004, 4.029008370204711]
ramp_means_carbon=[0.23192913053278295, 0.07597303082211063,133.999986524228, 22.68302, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 2.8654939556021e-10, 8.940621635999642,0.7999999909807196,4.749140535109004, 4.029008370204711]
#noramp_fit.optim_list=['E_0', 'k_0', 'Ru','Cdl', 'gamma','omega', 'phase','alpha']d
carbon_means=[0.2384234383845605, 2.8016906117813805, 14.009822731553978, 3.325500615749805e-05, 0.12980099583898141, -0.005247496549209238, 1.320027160165847e-10, 8.940856411874309, 5.513246206729991, 4.319227044942643, 0.5815664903172559]
noramp_fit.optim_list=['E_0','k_0', 'Ru',"Cdl","CdlE1", "CdlE2",'gamma', 'omega', "phase", "cap_phase", "alpha"]
test=noramp_fit.test_vals(carbon_means_2, "timeseries", test=False)
plt.plot(time_results,np.flip(test))
plt.plot(time_results, current_results)
plt.show()
desired_params=['E_0','k_0', 'Ru',"Cdl","CdlE1", "CdlE2",'gamma', "phase", "cap_phase", "alpha"]
plot_params=noramp_fit.pick_paramaters(carbon_means_2, desired_params)
noramp_fit.param_scanner(plot_params, desired_params, 0.2)
ramp_free_black=(noramp_fit.test_vals(carbon_means_2, "timeseries", test=False))
ramp_free_harmonics=harm_class.generate_harmonics(time_results, ramp_free_black)
#exp_harmonics_2=harm_class.generate_harmonics(time_results, test_2)
noramp_fit.optim_list=['E0_mean', "E0_std",'k_0', 'Ru',"Cdl", "CdlE1", "CdlE2",'gamma', 'omega', "alpha", "phase", "cap_phase"]

ramp_time=noramp_fit.test_vals(ramp_means_carbon, "timeseries", test=False)
ramp_harmonics=harm_class.generate_harmonics(time_results, ramp_time)
#harm_class.plot_harmonics(time_results, method="phased", Experimental=data_harmonics, Ramp_free=exp_harmonics, Ramped=ramp_harmonics)
harm_class.harmonics_and_voltages(noramp_fit.t_nondim(time_results),noramp_fit.e_nondim(voltage_results), folder, "phased", \
                            Experimental_harmonics=noramp_fit.i_nondim(data_harmonics), Sinusoidal_harmonics=noramp_fit.i_nondim(ramp_free_harmonics),\
                             Experimental_time_series=noramp_fit.i_nondim(current_results),Sinusoidal_time_series=noramp_fit.i_nondim(ramp_free_black))

dummy_times=np.linspace(0, 1, len(likelihood_func))
#noramp_fit.optim_list=['Ru', 'omega']
noramp_fit.dim_dict["alpha"]=0.5
noramp_fit.optim_list=['E_0','k_0',"Ru","Cdl","CdlE1", "CdlE2",'gamma', 'omega', "phase", "cap_phase"]
param_boundaries=np.zeros((2, noramp_fit.n_parameters()))
for i in range(0, noramp_fit.n_parameters()):
    param_boundaries[0][i]=param_bounds[noramp_fit.optim_list[i]][0]
    param_boundaries[1][i]=param_bounds[noramp_fit.optim_list[i]][1]
noramp_fit.define_boundaries(param_boundaries)
fourier_arg=likelihood_func
true_data=current_results
#if simulation_options["experimental_fitting"]==False:
#elif simulation_options["experimental_fitting"]==True:
    #fourier_arg=likelihood_func
    #true_data=current_results
noramp_fit.pass_extra_data(true_data, fourier_arg)
if simulation_options["likelihood"]=="timeseries":
    cmaes_problem=pints.SingleOutputProblem(noramp_fit, time_results, true_data)
elif simulation_options["likelihood"]=="fourier":
    cmaes_problem=pints.SingleOutputProblem(noramp_fit, dummy_times, fourier_arg)

score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(noramp_fit.optim_list))], [np.ones(len(noramp_fit.optim_list))])
noramp_fit.simulation_options["label"]="cmaes"
x0=abs(np.random.rand(noramp_fit.n_parameters()))#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
print len(x0), noramp_fit.n_parameters()
num_runs=20
score_mat=np.zeros(num_runs)
param_mat=np.zeros((len(noramp_fit.optim_list), num_runs))
for i in range(0, num_runs):
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
    print found_value
    cmaes_results=noramp_fit.change_norm_group(found_parameters, "un_norm")
    print list(cmaes_results)
    cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=True)
    plt.subplot(4,5,i+1)
    plt.plot(voltage_results, current_results)
    plt.plot(voltage_results, cmaes_time)
plt.show()


print list(cmaes_results)
#noramp_fit.simulate(found_parameters,time_results, normalise=True, likelihood="fourier", test=True )
cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=True)
cmaes_fourier=noramp_fit.test_vals(cmaes_results, likelihood="fourier", test=True)
print folder
print Method
plt.plot(time_results, true_data)
plt.plot(time_results, cmaes_time)
plt.show()
plt.plot(np.fft.ifft(cmaes_fourier))
plt.plot(np.fft.ifft(likelihood_func))
plt.show()

cmaes_harmonics=harm_class.generate_harmonics(time_results, cmaes_time)
true_vals=copy.deepcopy(cmaes_results)
for i in range(0, len(noramp_fit.optim_list)):
    plt.subplot(2,6, i+1)
    plt.plot(voltage_results, true_data)
    true_value=cmaes_results[i]
    for j in np.arange(0.5, 2.0, 0.5):
        cmaes_results[i]=true_value*j
        time_val=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=True)
        plt.plot(voltage_results, time_val, label=noramp_fit.optim_list[i]+"="+str(true_value*j))
    plt.legend()
    cmaes_results=true_vals
plt.show()
#hann=np.hanning(len(cmaes_time))
#f=np.fft.fftfreq(len(time_results), time_results[1]-time_results[0])
#Y1=np.multiply(hann, cmaes_time)
#y2=np.multiply(hann, current_results)
#Y1=np.fft.fft(Y1)
#y2=np.fft.fft(y2)
#exp_harmonics=harm_class.generate_harmonics(time_results, cmaes_time)
#harm_class.plot_harmonics(time_results, exp_harmonics, data_harmonics)
#plt.plot(f,abs(Y1))
#plt.plot(f,abs(y2))
#plt.show()

harm_class.plot_harmonics(time_results, cmaes_harmonics, data_harmonics,"abs", "numerical", "data")
cmaes_prediction=noramp_fit.simulate(found_parameters,frequencies, normalise=True, likelihood="timeseries", test=True )
print cmaes_results

#error=np.std(np.subtract(cmaes_prediction, likelihood_func))
"""
#error=np.std(np.subtract(cmaes_time, current_results))
mcmc_problem=pints.SingleOutputProblem(noramp_fit, dummy_times, likelihood_func)
#mcmc_problem=pints.SingleOutputProblem(noramp_fit, time_results, current_results)
updated_lb=np.append(cmaes_results*0.5, [0])#found_parameters[3]*0.97,
updated_ub=np.append(cmaes_results*1.5, [max(likelihood_func)])#found_parameters[3]*1.03,
updated_ub[1]=1e6
updated_boundaries=[updated_lb, updated_ub]
updated_boundaries=np.sort(updated_boundaries, 0)
noramp_fit.define_boundaries(updated_boundaries)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                updated_boundaries[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(cmaes_results, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
noramp_fit.label="MCMC"
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()
#pints.plot.trace(chains)
#plt.show()
means=[np.mean(chains[:, 5000:, x]) for x in np.arange(0, len(noramp_fit.optim_list))]
print means
print str(folder)
#mcmc_times=noramp_fit.simulate(means, time_results, "no", "timeseries", "yes")
#mcmc_harmonics=harm_class.generate_harmonics(time_results, mcmc_times)
#harm_class.plot_harmonics(time_results, mcmc_harmonics, data_harmonics,"abs")
#harm_class.plot_harmonics(time_results, mcmc_harmonics, data_harmonics, "negs")
#open_query=raw_input("save?")

filename=str(desired_length)+"_"+str(dec_amount)+"_"+str(lcv_3)+"."+folder[1:]
f=open(filename, "w")
np.save(f, chains)
f.close()
"""
#means=[2.85905314e-01, 5.86252081e+00, 1.19877032e-10, 4.19903721e-05, 4.78258908e+02, 8.94055732e+00]
#means2=[9.30117566e+03,1.16018837e-10,1.05469193e-05,1.63103403e+03,8.94085502e+00]
means=[2.36958426e-01, 9.99999958e+03, 1.58653602e+02, 1.50881707e-05, 6.37371838e-11, 8.94639042e+00,4.71238898038469+math.pi]#9.58653602e+02#8.56987286e+02#1.62644682e+02
means2=[2.36958426e-01, 9.99999958e+03, 1.58653602e+02, 1.50881707e-05, 6.37371838e-11, 8.94639042e+00,4.71238898038469]
inv_cdl_means_black=[0.022625846615525375, 282.2528765716474, 1565.32008743633, 7.233724535557822e-06, 0.0432328629721388, 0.00039664847059295294, 1.1720248951887567e-10, 8.940733606629868, 1.4114762214007537, 1.6792883395043496e-06]
inv_cdl_means_red=[0.22575364212262392, 197.9059067794697, 1939.9502453947625, 1.289377946400102e-05, 0.03672445932151236, -0.0003299037569737262, -1.1304027040442977e-05, 8.677149187960197e-11, 8.940937039777493, 1.5224066122900832, 1.3243217094538047e-08]
inv_cdl_means_plain=[0.22344327433678945, 254.09569103628488, 2984.336479477873, 5.488953773286254e-06, 0.030523464663436695, 0.0002271156823834275, 5.924975465640491e-11, 8.940811663124784, 1.4491256345114145, 4.729976917695294e-09]
ninv_cdl_means_red=[0.246621930911074, 78.45894300594209, 2818.6004999485217, 7.053324170917987e-07, -0.673031618742151, 0.09868327578110173, 1.4176466564885793e-10, 8.941170291777855, 2.034070948215475, 0.10000000056123999]
ninv_cdl_means_black=[0.3377987278522049, 0.6859606100145135, 1.8875103483928486e-07, 3.014005878620766e-17, -0.46299465435793596, 0.09248778871530525, 1.9193462079072354e-10, 8.940774114010635, 2.6487798785179644, 0.3622924808900895]
ninv_cdl_means_plain=[0.2660521816274963, 92.25120877213567, 2565.953228350421, 2.0728666764722445e-06, 0.09437714769471839, 0.0029509862465845887, 6.681522712570585e-11, 8.940948180628803, 1.8838080498194736, 0.10000000581935341]
means=[0.2, 100, 50, 1e-5, 0,0,1e-10, 8.94, 3*math.pi/2, 0.5]
GoldLarge1=[0.2634471594256482, 182.46968143921333, 251.52483060070875, 9.948258528983337e-05, 0.021593193177659842, 6.486044409236416e-10, 8.941736782729983, 1.6346028726977888, 0.10000000000199746]
GoldLarge2=[0.26039781003762164, 166.631058358149, 567.8337430471373, 4.3204616818471545e-05, 0.024366813205182414, 3.1255819860958073e-10, 8.941659017734578, 1.6575003335360676, 0.10000000000002386]
GoldLarge3=[0.2584206501377829, 157.11505190360802, 2749.660863756472, 8.802056383466946e-06, 0.026489956471885456, 6.70238627793266e-11, 8.94154301432433, 1.6692873967481876, 0.10000000058720873]
GC4_2=[0.2568519753260954, 384.6885351368988, 724.4994785328084, 2.9476802894874582e-05, 0.022962092313146165, 1.0351489419094205e-10, 8.940907182492527, 1.4050637398962589, 0.10000000000006122]
GC4_1=[0.257740555023939, 357.3609131669447, 382.2055036116924, 5.440677575193328e-05, 0.026386317796860403, 1.9319598326712751e-10, 8.94098189688349, 1.406746302896052, 0.1000000005558953]
#GC4_1=[0.2577409256040222, 357.41199962004214, 208.34836991645813, 9.981079168966267e-05, 0.02638451283800336, 3.544175619714202e-10, 8.940982084935797, 1.4067453971726527, 0.10000000009569851]
GC4_3=[0.25636968524906856, 393.3143126102424, 221.27196795533328, 9.809762640551237e-05, 0.021975801504555026, 3.453653918573443e-10, 8.940903530763926, 1.4078216880545553, 0.10000000000000066]
means=[0.22, 500, 200, 5e-5, 1e-10, 8.94, 3*math.pi/2, 0.5]
means=[0.217740555023939, 357.3609131669447, 382.2055036116924, 5.440677575193328e-05, 0, 1.9319598326712751e-10, 8.94098189688349, 3/2*math.pi, 0.9000000005558953]
ramp_means=[0.2112423897575218, 0.07357984171424121, 126.17588866800264, 24.24736830211999, 9.999999999876686e-07, 0.11886527586092166, 0,1.3175840584818795e-10, 8.940448789535177, 0.599999999999731,4.749140535109004, 4.749140535109004]
ramp_means=[0.21200070197828386, 0.06888959981526988, 133.96563492653507, 40.08177557223102, 3.226207450320691e-06, -0.021487125154184827, 0.0017931237883237632, 1.2669148148700962e-10, 8.940448789535177, 0.7999999999999661,4.749140535109004, 4.028465609498452]
ramp_means_carbon_1=[0.24192913053278295, 9999.999986524228, 1.0428644292961376e-08, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 3.5654939556021e-10, 8.940621635999642, 4.749140535109004, 4.029008370204711]
ramp_means_carbon=[0.23192913053278295, 0.07597303082211063,133.999986524228, 22.68302, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 2.8654939556021e-10, 8.940621635999642,0.7999999909807196,4.749140535109004, 4.029008370204711]
#best_so_far="[4.57076403e-01 2.76438997e-02 1.00989565e-01 4.83961049e-06 1.43033271e-01]"
#best_so_far_red=[0.2876, 10000, 660, 0.000094, 1.47e-10]
#best_high_freq=[0.237, 13.5, 120, 0.0020, 4.1e-10]
noise_val=0.02
noise_max=max(test)*noise_val
noise=np.random.normal(0,noise_max, len(test))
synthetic_data=np.add(test, noise)
fourier_test1=noramp_fit.kaiser_filter(synthetic_data)
test_data=np.fft.ifft(likelihood_func)
fourier_test=np.fft.ifft(fourier_test1)
L=len(test)
time_len=range(0, len(time_results))
f_len=np.linspace(0, time_results[-1], len(fourier_test))
time_plot=np.interp(f_len, time_len, time_results)
hann=np.hanning(L)
