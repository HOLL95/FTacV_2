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

dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Carbon"#"/Gold/Large"
Method ="GC4_1"#GoldLarge_1"
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

desired_length=int(2e4)
dec_amount=8
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
voltage_results=results2[0::dec_amount, 1]
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.09, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-5, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1e-10,          # (surface coverage per unit area)
    'k_0': 10, #(reaction rate s-1)
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
    "time_end": None,
    'num_peaks': 50
}
solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
likelihood_options=["timeseries", "fourier"]
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "test":False,
    "likelihood":likelihood_options[1],
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
    "bounds_val":2000,
    "signal_length":int(1e4),
}
param_list['E_0']=2.36504746e-01#(param_list['E_reverse']-param_list['E_start'])/2
noramp_fit=single_electron(param_list, simulation_options, other_values)
time_results=noramp_fit.other_values["experiment_time"]
current_results=noramp_fit.other_values["experiment_current"]
voltage_results=noramp_fit.other_values["experiment_voltage"]
likelihood_func=noramp_fit.kaiser_filter(current_results)
test_voltages=noramp_fit.define_voltages()
noramp_fit.optim_list=[]
noramp_fit.bounds_val=max(current_results)/10
frequencies=noramp_fit.frequencies
noramp_fit.pass_extra_data(current_results, likelihood_func)
test1=noramp_fit.simulate([],frequencies, normalise=False, likelihood="timeseries", test=False)#
for lcv_1 in np.arange(0.05, 0.15, 0.01):
    for lcv_2 in range(0,3):
        dec_amount=8
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
            'area': lcv_1, #(electrode surface area cm^2)
            'Ru': 1.0,  #     (uncompensated resistance ohms)
            'Cdl': 1e-5, #(capacitance parameters)
            'CdlE1': 0,#0.000653657774506,
            'CdlE2': 0,#0.000245772700637,
            'CdlE3': 0,#1.10053945995e-06,
            'gamma': 1e-10,          # (surface coverage per unit area)
            'k_0': 10, #(reaction rate s-1)
            'alpha': 0.5,
            'sampling_freq' : (1.0/200),
            'phase' : 3*(math.pi/2),
            "time_end": None,
            'num_peaks': 50
        }
        solver_list=["Bisect", "Brent minimisation", "Newton-Raphson", "inverted"]
        likelihood_options=["timeseries", "fourier"]
        simulation_options={
            "no_transient":False,
            "numerical_debugging": False,
            "experimental_fitting":True,
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
            "bounds_val":2000,
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
        noramp_fit.bounds_val=max(current_results)/10
        frequencies=noramp_fit.frequencies
        noramp_fit.pass_extra_data(current_results, likelihood_func)
        test1=noramp_fit.simulate([],frequencies, normalise=False, likelihood="timeseries", test=False)#

        param_bounds={
            'E_0':[0.1, 0.4],#[param_list['E_start'],param_list['E_reverse']],
            'omega':[0.98*param_list['omega'],1.02*param_list['omega']],#8.88480830076,  #    (frequency Hz)
            'Ru': [0, 3e3],  #     (uncompensated resistance ohms)
            'Cdl': [0,1e-4], #(capacitance parameters)
            'CdlE1': [-2,2],#0.000653657774506,
            'CdlE2': [-0.1,0.1],#0.000245772700637,
            'CdlE3': [-0.05,0.05],#1.10053945995e-06,
            'gamma': [1e-11,1e-9],
            'k_0': [0, 1e4], #(reaction rate s-1)
            'alpha': [0.1, 0.9],
            'phase' : [0, 2*math.pi]
        }
        noramp_fit.optim_list=['E_0', 'k_0', 'Ru','Cdl', 'CdlE1','CdlE2','gamma','omega', 'phase']
        #'E_0', 'k_0' 'Ru', 'Cdl','gamma'
        harm_class=harmonics(other_values["harmonic_range"], noramp_fit.nd_param.omega*noramp_fit.nd_param.c_T0, 0.1)
        data_harmonics=harm_class.generate_harmonics(time_results, current_results)
        #means=[2.85905314e-01, 5.86252081e+00, 1.19877032e-10, 4.19903721e-05, 4.78258908e+02, 8.94055732e+00]
        #means2=[9.30117566e+03,1.16018837e-10,1.05469193e-05,1.63103403e+03,8.94085502e+00]
        means=[param_list[key] for key in noramp_fit.optim_list]
        means=[2.36958426e-01, 9.99999958e+03, 1.58653602e+02, 1.50881707e-05, 6.37371838e-11, 8.94639042e+00,4.71238898038469+math.pi]#9.58653602e+02#8.56987286e+02#1.62644682e+02
        means2=[2.36958426e-01, 9.99999958e+03, 1.58653602e+02, 1.50881707e-05, 6.37371838e-11, 8.94639042e+00,4.71238898038469]
        inv_cdl_means_black=[0.022625846615525375, 282.2528765716474, 1565.32008743633, 7.233724535557822e-06, 0.0432328629721388, 0.00039664847059295294, 1.1720248951887567e-10, 8.940733606629868, 1.4114762214007537, 1.6792883395043496e-06]
        inv_cdl_means_red=[0.22575364212262392, 197.9059067794697, 1939.9502453947625, 1.289377946400102e-05, 0.03672445932151236, -0.0003299037569737262, -1.1304027040442977e-05, 8.677149187960197e-11, 8.940937039777493, 1.5224066122900832, 1.3243217094538047e-08]
        inv_cdl_means_plain=[0.22344327433678945, 254.09569103628488, 2984.336479477873, 5.488953773286254e-06, 0.030523464663436695, 0.0002271156823834275, 5.924975465640491e-11, 8.940811663124784, 1.4491256345114145, 4.729976917695294e-09]
        ninv_cdl_means_red=[0.246621930911074, 78.45894300594209, 2818.6004999485217, 7.053324170917987e-07, -0.673031618742151, 0.09868327578110173, 1.4176466564885793e-10, 8.941170291777855, 2.034070948215475, 0.10000000056123999]
        ninv_cdl_means_black=[0.3377987278522049, 0.6859606100145135, 1.8875103483928486e-07, 3.014005878620766e-17, -0.46299465435793596, 0.09248778871530525, 1.9193462079072354e-10, 8.940774114010635, 2.6487798785179644, 0.3622924808900895]
        ninv_cdl_means_plain=[0.2660521816274963, 92.25120877213567, 2565.953228350421, 2.0728666764722445e-06, 0.09437714769471839, 0.0029509862465845887, 6.681522712570585e-11, 8.940948180628803, 1.8838080498194736, 0.10000000581935341]
        means=[0.2, 100, 50, 1e-5, 0,0,1e-10, 8.94, 3*math.pi/2, 0.5]
        #means=[3.47301042e-01, 2.24569221e+03, 4.50000000e+05, 1.00000000e-06, 7.70422832e-10, 8.92022169e+00]
        #means=means4_8
        alpha_range=[0.05, 0.2,0.5, 0.7, 0.95]
        noramp_fit.simulation_options["likelihood"]="timeseries"
        noramp_fit.optim_list=['E_0', 'k_0', 'Ru','Cdl', 'CdlE1','CdlE2','gamma','omega', 'phase','alpha']
        test2=noramp_fit.simulate(inv_cdl_means_black,frequencies, normalise=False, likelihood="timeseries", test=False )
        noramp_fit.optim_list=['E_0', 'k_0', 'Ru','Cdl', 'CdlE1','CdlE2','gamma','omega', 'phase','alpha']
        test=noramp_fit.simulate(means,frequencies, normalise=False, likelihood="timeseries", test=False )
        #noramp_fit.simulation_options["likelihood"]="fourier"
        #test=np.flip(np.array(test)*-1)
        #test2=noramp_fit.simulate(means2,frequencies, "no", "timeseries", "no" )
        noise_val=0.02
        noise_max=max(test)*noise_val
        noise=np.random.normal(0,noise_max, len(test))
        synthetic_data=np.add(test, noise)
        #plt.plot(synthetic_data)
        #plt.show()
        #plt.plot(likelihood_func)
        #plt.show()
        #fourier_test1=noramp_fit.simulate(means,frequencies, "no", "fourier", "no" )
        fourier_test1=noramp_fit.kaiser_filter(synthetic_data)
        test_data=np.fft.ifft(likelihood_func)
        fourier_test=np.fft.ifft(fourier_test1)
        L=len(test)
        time_len=range(0, len(time_results))
        f_len=np.linspace(0, time_results[-1], len(fourier_test))
        time_plot=np.interp(f_len, time_len, time_results)
        hann=np.hanning(L)
        #plt.plot(f_len, (fourier_test), label="numerical")
        #plt.plot(f_len, (test_data), label="data", linestyle="-", linewidth=3, alpha=0.7)
        #plt.legend()
        #plt.show()
        #noramp_fit.simulate(means,frequencies, "no", "fourier", "yes" )
        exp_harmonics=harm_class.generate_harmonics(time_results, test)

        #harm_class.plot_harmonics(time_results, exp_harmonics,data_harmonics, "phased")
        #harm_class.harmonics_and_time(time_results, exp_harmonics, np.flip(test), current_results, data_harmonics, "numerical", "data")
        #exp_harmonics=harm_class.generate_harmonics(time_results, test2)
        #harm_class.harmonics_and_time(time_results, exp_harmonics, test2, current_results, data_harmonics, "numerical", "data")

        dummy_times=np.linspace(0, 1, len(likelihood_func))
        #noramp_fit.optim_list=['Ru', 'omega']
        noramp_fit.optim_list=['E_0','k_0', 'Ru', 'Cdl','CdlE1','CdlE2','gamma', 'omega', 'phase', "alpha"]
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
        for i in range(0, 1):
            found_parameters, found_value=pints.optimise(
                                                        score,
                                                        x0,
                                                        boundaries=CMAES_boundaries,
                                                        method=pints.CMAES
                                                        )
        cmaes_results=noramp_fit.change_norm_group(found_parameters, "un_norm")
        print list(cmaes_results)
        noramp_fit.simulate(found_parameters,time_results, normalise=True, likelihood="fourier", test=True )
        noramp_fit.simulation_options["likelihood"]="timeseries"
        cmaes_time=noramp_fit.simulate(found_parameters,time_results, normalise=True, likelihood="timeseries", test=False )
        print folder
        plt.plot(time_results, true_data)
        plt.plot(time_results, cmaes_time)
        plt.show()
        plt.plot(voltage_results, true_data)
        plt.plot(voltage_results, cmaes_time)
        plt.show()
        plt.plot(voltage_results, true_data)
        plt.plot(voltage_results, np.flip(cmaes_time, axis=0))
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
        cmaes_harmonics=harm_class.generate_harmonics(time_results, cmaes_time)
        harm_class.plot_harmonics(time_results, cmaes_harmonics, data_harmonics,"phased", "numerical", "data")
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
        #best_so_far="[4.57076403e-01 2.76438997e-02 1.00989565e-01 4.83961049e-06 1.43033271e-01]"
        #best_so_far_red=[0.2876, 10000, 660, 0.000094, 1.47e-10]
        #best_high_freq=[0.237, 13.5, 120, 0.0020, 4.1e-10]
