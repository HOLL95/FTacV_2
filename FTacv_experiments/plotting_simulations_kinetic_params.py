import numpy as np
import matplotlib.pyplot as plt
import math
import time
start=time.time()
from single_e_class_unified import single_electron
from ramped_setup import FTACV_initialisation
from harmonics_plotter import harmonics
print "import", time.time()-start
types=["current", "voltage"]
start=time.time()
noramp_startup=FTACV_initialisation(experimental_fitting=True, file_dict={"GC4_1_cv":types, "GC4_2_cv":types, "GC4_3_cv":types}, dec_amount=4)
ramp_startup=FTACV_initialisation(experimental_fitting=True, file_dict={"GC4_1_ramp_cv":types,}, dec_amount=64)
print ramp_startup.current_results
print "read", time.time()-start
simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":40,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"sinusoidal",
    "label": "MCMC",
    "optim_list":[]
}
ramped_simulation_options={
    "no_transient":False,
    "numerical_debugging": False,
    "experimental_fitting":True,
    "dispersion":False,
    "dispersion_bins":40,
    "test": False,
    "likelihood":"timeseries",
    "numerical_method": "Brent minimisation",
    "method":"ramped",
    "label": "MCMC",
    "optim_list":[]
}
noramp_other_values={
    "filter_val": 0.5,
    "harmonic_range":range(3,7,1),
    "experiment_time": noramp_startup.time_results["GC4_1_cv"],
    "experiment_current":noramp_startup.current_results["GC4_1_cv"],
    "experiment_voltage":noramp_startup.voltage_results["GC4_1_cv"],
    "bounds_val":20,
    "signal_length":int(2e4),
}

ramped_other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,9,1),
    "experiment_time": ramp_startup.time_results["GC4_1_ramp_cv"],
    "experiment_current": ramp_startup.current_results["GC4_1_ramp_cv"],
    "experiment_voltage":ramp_startup.voltage_results["GC4_1_ramp_cv"],
    "bounds_val":20,
    "signal_length":int(5e5),
}

noramp_simulations=single_electron(noramp_startup.generic_noramp_params, simulation_options, noramp_other_values)
ramped_simulations=single_electron(ramp_startup.generic_ramped_params, ramped_simulation_options, ramped_other_values)
harm_class=harmonics(noramp_other_values["harmonic_range"],noramp_simulations.nd_param.omega, 0.2)
ramped_simulations.def_optim_list(["k_0", "omega"])
omegas=[math.pi*2**x for x in range(1, 6)]
k_s=np.flip([10**x for x in range(1, 5)])
fig, ax =plt.subplots(1, len(omegas))
for i in range(0, len(omegas)):
    for j in range(0, len(k_s)):
        ax[i].plot(ramped_simulations.time_vec, ramped_simulations.i_nondim(1e3*ramped_simulations.test_vals([k_s[j], omegas[i]], "timeseries")), label=k_s[j])
    ax[i].set_title('$\\omega$='+ str(omegas[i]/math.pi)+"$\pi$")
    ax[i].set_xlabel("Time(s)")
    ax[i].set_ylabel("Current(mA)")
    ax[i].legend()
plt.show()
