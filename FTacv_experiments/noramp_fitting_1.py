import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
from single_e_class_noramp  import single_electron
from harmonics_plotter import harmonics
from matplotlib.widgets import Slider, Button, RadioButtons
import pints
import pints.plot
import copy
import time
params_for_opt=[]
def file_opener(files, file_path, file_dict, dec_amount):
    voltage_results={}
    time_results={}
    current_results={}
    experiments=file_dict.keys()
    for data in files:
        for keys in experiments:
            Method=keys
            type=file_dict[keys]
            if (Method in data)  and (type[0] in data):
                file=path+"/"+data
                results=np.loadtxt(file)
                current_results[Method]=results[0::dec_amount, 1]
                time_results[Method]=results[0::dec_amount, 0]#
                break
            elif (Method in data)  and (type[1] in data):
                file=path+"/"+data
                results=np.loadtxt(file)
                voltage_results[Method]=results[0::dec_amount, 1]#
                break
    return voltage_results, current_results, time_results
folder_options=["/Black", "/Red", "/Plain", "/Carbon", "/Gold/Large", "/Gold/Small"]
dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Carbon"#"/Black"#"/Gold/Large"#
Method ="PostSinusoidal"
type="asA_4"
path=dir_path+data_path+folder
files= os.listdir(path)
types=["current", "voltage"]
file_dict={"GC4_1_cv":types, "GC4_2_cv":types, "GC4_3_cv":types}
voltages, currents, times=file_opener(files, path, file_dict, 4)
filenames=sorted(file_dict.keys())
for files in filenames:
    plt.plot(voltages[files], currents[files], label=files)
plt.legend()
plt.show()
time_results=times["GC4_1_cv"]
voltage_results=voltages["GC4_1_cv"]
#plt.plot(time_results, current_results)
#plt.show()

de=300e-3
estart=260e-3-de
ereverse=estart+2*de
file="GC4_1_cv"
param_list={
    "E_0":0.2,
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.07, #(electrode surface area cm^2)
    'Ru': 1.0,  #     (uncompensated resistance ohms)
    'Cdl': 1e-6, #(capacitance parameters)
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
    "dispersion_bins":20,
    "test": False,
    "likelihood":likelihood_options[0],
    "numerical_method": solver_list[1],
    "label": "MCMC",
    "optim_list":[]
}
other_values={
    "filter_val": 0.5,
    "harmonic_range":range(1,9,1),
    "experiment_time": times[file],
    "experiment_current": currents[file],
    "experiment_voltage":voltages[file],
    "bounds_val":20,
    "signal_length":int(2e4),
}
param_bounds={
    'E_0':[0.2, 0.3],#[param_list['E_start'],param_list['E_reverse']],
    'omega':[0.95*param_list['omega'],1.05*param_list['omega']],#8.88480830076,  #    (frequency Hz)
    'Ru': [0, 500],  #     (uncompensated resistance ohms)
    'Cdl': [0,1e-4], #(capacitance parameters)
    'CdlE1': [-0.05,0.15],#0.000653657774506,
    'CdlE2': [-0.008,0.008],#0.000245772700637,
    'CdlE3': [-0.01,0.01],#1.10053945995e-06,
    'gamma': [1e-11,1e-9],
    'k_0': [1, 500], #(reaction rate s-1)
    'alpha': [0.4, 0.6],
    "cap_phase":[0, 2*math.pi],
    "E0_mean":[0.15, 0.3],
    "E0_std": [0.01, 0.2],
    "k0_shape":[0,5],
    "k0_loc":[1, 1e3],
    "k0_scale":[0,2e3],
    "k0_range":[1e2, 1e4],
    'phase' : [0, 2*math.pi]
}
#(param_list['E_reverse']-param_list['E_start'])/2
noramp_fit=single_electron(param_list, simulation_options, other_values)
noramp_fit.define_boundaries(param_bounds)
time_results=noramp_fit.other_values["experiment_time"]
current_results=noramp_fit.other_values["experiment_current"]
voltage_results=noramp_fit.other_values["experiment_voltage"]
for i in range(0, len(filenames)):
    voltages[filenames[i]]=voltages[filenames[i]][noramp_fit.time_idx:other_values["signal_length"]]
    currents[filenames[i]]=currents[filenames[i]][noramp_fit.time_idx:other_values["signal_length"]]
    times[filenames[i]]=times[filenames[i]][noramp_fit.time_idx:other_values["signal_length"]]
likelihood_func=noramp_fit.kaiser_filter(current_results)
test_voltages=noramp_fit.define_voltages()
frequencies=noramp_fit.frequencies
noramp_fit.pass_extra_data(current_results, likelihood_func)

unit_dict={
    "E_0": "V",
    'E_start': "V", #(starting dc voltage - V)1254684939476, 1.4485673957633017e-10, 8.94084836566341, 4.898877306283271, 4.379038086713167, 0.6999999953311751]

    'E_reverse': "V",
    'omega':"Hz",#8.88480830076,  #    (frequency Hz)
    'd_E': "V",   #(ac voltage amplitude - V) freq_range[j],#
    'v': r'$s^{-1}$',   #       (scan rate s^-1)
    'area': r'$cm^{2}$', #(electrode surface area cm^2)
    'Ru': r'$\Omega$',  #     (uncompensated resistance ohms)
    'Cdl': "F", #(capacitance parameters)
    'CdlE1': "",#0.000653657774506,
    'CdlE2': "",#0.000245772700637,
    'CdlE3': "",#1.10053945995e-06,
    'gamma': r'$mol$ $cm^{-2}$',
    'k_0': r'$s^{-1}$', #(reaction rate s-1)
    'alpha': "",
    "E0_mean":"V",
    "E0_std": "V",
    "k0_shape":"",
    "k0_loc":"",
    "k0_scale":"",
    "cap_phase":"rads",
    'phase' : "rads",
}

#'E_0', 'k_0' 'Ru', 'Cdl','gamma'
harm_class=harmonics(other_values["harmonic_range"], noramp_fit.nd_param.omega*noramp_fit.nd_param.c_T0, 0.1)
data_harmonics=harm_class.generate_harmonics(time_results, current_results)

#carbon_means=[0.24432865369586831, 0.03718591974278576, 5266.813117359024, 2.3351875006109577e-06, 3.904567925683219e-05, 0.020282592797868926,0, 2.886466778833776e-10, 8.940621635999642, 0.10000000291906254,]
carbon_means=[0.2629196040175071, 0.010000000426528552, 242.94186943883, 793.6555564492475, 4.020715996280745e-05, 0.012041277069243517,0, 1.702054870510513e-10, 8.94103486982429, 0.8999999976170034, 5.132595831851037, 4.446662466382718]
black_means=[0.20504878429957712, 0.04692985835905884, 773.0039074468887, 1.1172494386860095e-06, 6.253912022948444e-06, 0.46590284463560927, -0.020672008906663236, 9.665438282298918e-11, 8.94055300929529, 0.10000008017506497, 4.749140535109004, 4.029008370204711]
black_means_2=[0.22026089333976873, 0.04776183826475387, 1226.0003897156193, 2.6091841820962103e-10, 6.311657164346574e-06, 0.4861839637949368, -0.021712454485320096, 9.470125832724226e-11, 8.940448789535177, 0.10000000178756306, 4.742083569735882, 4.028465609498452]
black_means_3=[0.19813523414292655, 1.4877187307507795, 4.4556439040686615e-08, 1.1663407409020574e-06, 1.775725496685462, -0.07670261113515547, 1.2372296727209742e-10, 8.940758763885253, 5.622404458746398, 4.172942575228049, 0.5999999999990152]
carbon_means_2=[0.24437320720249422, 0.9641045582560264, 4.4556439040686615e-08, 3.3222199217955305e-05, 0.09074940228372252, -0.0034091889329414538, 1.6095034096343578e-10, 8.940780387771508, 5.698538006754413, 4.317679292947246, 0.5417098085713957]
carbon_means_3_alpha=[0.26689122715624614, 97.72540069109925, 601.2211952415878, 3.379513907637756e-05, 0.10070581762902808, -0.0038788288927744663, 1.4557809348768636e-10, 8.940779258804525, 5.149017610779459, 4.450212193055553, 0.5]
carbon_means_4=[0.26380081065342453, 71.68553415203698, 473.6569599754532, 3.2962209403627705e-05, 0.11560665350196508, -0.00447690722058372, 1.4277760997483275e-10, 8.940785781966794, 5.139463149738304, 4.428851337096966, 0.5999999999355239]
carbon_means_5_1000_k=[0.26764402303397383, 761.2191987917648, 3.450250003362729e-05, 0.08536733387749873, -0.003128578066079607, 1.5151605215241677e-10, 8.94061684149824, 5.100339061687143, 4.488467382773432, 0.5999999999903135]
carbon_means_6_20_ru=[ 2.45827239e-01, 1.038,20,  3.31526662e-05,  9.25625459e-02, -3.49813855e-03,  1.59975902e-10,  8.94073365e+00,  5.69475975e+00, 4.32091602e+00,  5.39205070e-01]
carbon_means_abs=[0.25839275077930485, 20.108646168645084, 570.1495572183396,3.2962209403627705e-05, 0.11560665350196508, -0.00447690722058372, 1.425656934986536e-10,8.940620380889355, 0.41318249114045236, 3*math.pi/2, 3*math.pi/2]
carbon_means_abs_all_peaks=[0.2712720627064147, 53.00972798907347, 126.43581153952566, 3.160962997490096e-05, 0.14647947194537103, -0.005848319334033306, 1.0072445202882476e-10, 8.940709552792356, 0.5488793068724522]
carbon_means_abs_all_peaks=[0.2528081478002836, 217.56996476308117, 529.7297882492157, 3.0166883670561975e-05, 0.14790906812787918, -0.005389008883307372, 1.1330727195554263e-10, 8.940664568706621, 0.40000000111206985]
another_constrained_fit=[0.25007589734382935, 0.04947794128755774, 20.83219016081861, 0.0002356574033788102, 3.361058693410845e-06, 1.9999995197072655, -0.09693566779590666, 3.665491921443447e-10, 8.941522982029218, 0.49515486096231986, 3.8227156574615946, 5.2169580962884226]
carbon_means_disped=[0.23542237635148117, 0.004406518077359679, 1.6338863524591192, 80.0000000070073, 3.24102216543367e-05, 0.10683350182222259, -0.0041305496741190365, 1.531049485509304e-10, 8.94090414236465, 5.635236560276772, 4.327586365963918, 0.5420683576660721]
carbon_means_disped_gc41a=[0.2507124585192858, 0.012810161448688611, 81.25609725526283,100, 3.4455963089647964e-05, 0.08570032724521258, -0.0026762829980299324, 1.4591816196158548e-10, 8.940730783653084, 4.965641139765749, 4.375333040247454, 0.7811192097028208]
carbon_means_disped_gc42a=[0.25150277022331746, 0.014115779656223786, 78.20813556543256,58.804257252122355, 3.501508623604994e-05, 0.08815081354697163, -0.0030182380193618527, 1.3910336262835592e-10, 8.941274817456382, 4.946748053531174, 4.3670012650213295, 0.7663320409205243]
carbon_means_disped_gc43a=[0.2508929602137271, 0.014962342550521588,  74.5650170269896,13.194027307777567, 3.51489903928365e-05, 0.0898311438945581, -0.0031909694983879135, 1.3515246782381657e-10, 8.941160154930586, 4.933768379207287, 4.358597908430476, 0.7621639225400778]
carbon_means_disped_gc42=[0.25983027085247, 78.20813556543256, 100.00001599685609, 3.501508623604994e-05, 0.08815081354697163, -0.0030182380193618527, 1.3910336262835592e-10, 8.941274817456382, 4.946748053531174, 4.3670012650213295, 0.7785371130874895]
carbon_means_disped_gc43=[0.25947921946089975, 74.5650170269896, 100.0000011034549117, 3.51489903928365e-05, 0.0898311438945581, -0.0031909694983879135, 1.3515246782381657e-10, 8.941160154930586, 4.933768379207287, 4.358597908430476, 0.7761838774833711]
higher_sampling_freq=[0.25699207481665864, 0.05555428879890889, 414.56746021325034, 325.8031006821857, 3.519141869992617e-05, 0.07717238077393472, -0.002727432717977485, 1.4690610585539404e-10, 8.940943575263907, 4.896585013895377, 4.391652613525725, 0.8999999860779]
higher_sampling_freq_2=[0.24936776983318537, 0.06166238798441354, 65.95536643340846, 50.00000000744157, 3.403237947373622e-05, 0.09662067753703223, -0.003521905241046774, 1.3876360325400223e-10, 8.940969817645122, 4.896305529543139, 4.339663409102047, 0.5766195385556642]
carbon_means_disped_gc41a=[0.24673357615129932, 0.06431065804807556, 61.288851318074876, 1.0000002125554084, 3.427115851232685e-05, 0.08880338471929529, -0.0033551664287031473, 1.3976795665535952e-10, 8.940855717890603, 4.877449155749592, 4.3314073608533965, 0.572544218603636]
gc4_2=[0.24784915687346862, 0.06316612861983209, 58.030622008233514, 1.0000000612187943, 3.405156321210489e-05, 0.09670855824995263, -0.0035129837977470983, 1.3801359615989042e-10, 8.94093762162615, 4.882482708448655, 4.3317118352679245, 0.5746764358057844]
gc4_3=[0.24784850701480876, 0.06316545250773771, 58.03170289032332, 1.0000000371416478, 3.405151295634379e-05, 0.0967118100534305, -0.0035131450128175305, 1.380117209701837e-10, 8.940937368486582, 4.882485181424615, 4.33171345860357, 0.5746856092254882]
carbon_means_disped_gc42a=[0.25386735650141046, 0.05718314856725754, 220.85131424751998, 263.61923792350825, 3.4658312714307376e-05, 0.08305864277903698, -0.003101254684939476, 1.4485673957633017e-10, 8.94084836566341, 4.898877306283271, 4.379038086713167, 0.6999999953311751]
gc4_3_alternate=[0.24517090736930636, 0.056720112426872436, 223.12968601072467, 269.4152026256152, 3.467489815250535e-05, 0.08073629430570951, -0.003038015139181096, 1.4366133260790924e-10, 8.940854696401248, 4.9005338607773, 4.378474411191984, 0.6999999999960397]
gc4_1_alternate=[0.24707364221482686, 0.05611157233858666, 215.93921919728885, 261.0012935458973, 3.457664799584686e-05, 0.08731475943652235, -0.00315097559270702, 1.4305526845381745e-10, 8.940962054668862, 4.897135123391333, 4.379405055574097, 0.6999999999984596]
gc4_2_alternate=[0.24581388949144, 0.056437609561567396, 220.84863259145277, 268.249787363191, 3.465024308699503e-05, 0.08279171693559563, -0.003087944274776373, 1.4514161563358487e-10, 8.940879185090578, 4.900935582894143, 4.379408688250694, 0.6999999997489996]

gc4_2_low_ru=[0.2390424106953814, 0.0634126339128772, 59.56331759126466, 1.0000000052336186, 3.421333538191425e-05, 0.09142331002140631, -0.0034273912843145503, 1.3806680869872118e-10, 8.940855654906262, 4.881140108985157, 4.332212053411549, 0.5726231710329381]
gc4_1_low_ru=[0.23953731576330606, 0.06251048824509087, 54.80881131862505, 1.0, 3.401227622853829e-05, 0.09658746281397658, -0.0035025931533504253, 1.4039829708596723e-10, 8.940959682968415, 4.88941265694472, 4.331036700670979, 0.5725681192157199]
gc4_3_low_ru=[0.23953952896099746, 0.06432508696935929, 60.84131849364619, 1.0, 3.42647900201487e-05, 0.08893011038618533, -0.0033602097013514043, 1.3972124127562154e-10, 8.94084399186188, 4.8779396488202185, 4.331377454498156, 0.5724578682161762]

gc4_2_b=[0.2390418691618472, 0.06341251668841151, 59.55992622848135, 1.0000000065311607, 3.42132227450153e-05, 0.09142523880452204, -0.003427483907861838, 1.4006677642620296e-10, 8.940855668037237, 4.881147386014593, 4.33221437864721, 0.5726264177402156]

carbon_means_disped=carbon_means_disped_gc42
#carbon_means_disped=[0.26309322667095436, 0.08411919477796326, 999.9999999749896, 5.005163812028062e-10, 2.7249704365065652e-05, 0.0014108474638033497, 0.0009417539364441913, 2.1916294751016806e-10, 8.94090414236465, 3*math.pi/2, 3*math.pi/2,0.7999999999895195]
#carbon_means_disped=[0.24042307692307694, 0.05995192307692309, 81.25609725526283, 100.0, 3.0169689037255904e-05, 0.020629272483467714, -0.0002036780323021597, 1.4591816196158548e-10, 8.959036763806955, 0.0,0.0, 0.5990384615384616]
#carbon_means_disped_2=[0.252803269870875, 0.07188314305956368, 88.46116906534513, 149.00689828448213, 3.64063523054103e-06, -0.13778150113414067, -0.0050711994601424434, 1.9534787417005916e-10, 8.94090414236465, 3*math.pi/2, 3*math.pi/2, 0.2570421859167864]
#carbon_means_disped=[0.26309322667095436, 0.08411919477796326, 999.9999999749896, 5.005163812028062e-10, 2.7249704365065652e-05, 0.0014108474638033497, 0.0009417539364441913, 2.1916294751016806e-10, 8.940730783653084, 3*math.pi/2, 3*math.pi/2, 0.5]
#carbon_means_disped=[0.2546334543150689, 0.03323126819215362, 116.24098132960108,144.78687238133318,3.4422425270083556e-05, 0.09197677289665342, -0.002971108200729257, 1.4130483196036997e-10,8.94090414236465, 4.947246577563367, 4.361083156674927, 0.5]
#carbon_means_abs=[0.25007589734382935, 20.83219016081861, 0.0002356574033788102, 3.361058693410845e-06, 1.9999995197072655, -0.09693566779590666, 3.665491921443447e-10, 0.49515486096231986, 3.8227156574615946, 5.2169580962884226]
#0.2547817849499173, 0.03183289363921536,


random_params=[0.25, 0.05, 100,100, 3e-5, 0.1, -0.005, 1e-10, 3*math.pi/2, 0.5]
noramp_fit.dim_dict["phase"]=3*math.pi/2
noramp_fit.dim_dict["cap_phase"]=3*math.pi/2
desired_params=['E0_mean', "E0_std",'k_0', 'Ru',"Cdl", "CdlE1", "CdlE2",'gamma', "cap_phase","alpha"]
noramp_fit.def_optim_list(desired_params)
#noramp_fit.param_scanner(random_params, desired_params, unit_dict, 0.2, "", boundaries=True)
carbon_means_8_cdl_only=[2.62877094e-01, 9.99984499e+03, 1.15828858e-06, 1.88609404e-05,7.79829258e-10, 8.94036654e+00, 4.79419627e+00, 5.99999976e-01]
ramp_means=[0.21200070197828386, 0.06888959981526988, 133.96563492653507, 40.08177557223102, 3.226207450320691e-06, -0.021487125154184827, 0.0017931237883237632, 1.2669148148700962e-10, 8.940448789535177, 0.7999999999999661,4.749140535109004, 4.028465609498452]
ramp_means_carbon_1=[0.24192913053278295, 9999.999986524228, 1.0428644292961376e-08, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 3.5654939556021e-10, 8.940621635999642, 4.749140535109004, 4.029008370204711]
ramp_means_carbon=[0.23192913053278295, 0.07597303082211063,133.999986524228, 22.68302, 9.999999999710656e-06, -0.19135843729198476, 0.012883589352296436, 2.8654939556021e-10, 8.940621635999642,0.7999999909807196,4.749140535109004, 4.029008370204711]
#noramp_fit.optim_list=['E_0', 'k_0', 'Ru','Cdl', 'gamma','omega', 'phase','alpha']d

carbon_means=[0.2384234383845605, 2.8016906117813805, 14.009822731553978, 3.325500615749805e-05, 0.12980099583898141, -0.005247496549209238, 1.320027160165847e-10, 8.940856411874309, 5.513246206729991, 4.319227044942643, 0.5815664903172559]
noramp_fit.def_optim_list(['E0_mean', "E0_std",'k_0', 'Ru',"Cdl", "CdlE1", "CdlE2",'gamma',"omega", "phase","cap_phase","alpha",])

noramp_fit.simulation_options["dispersion_bins"]=20
gc4_1=(noramp_fit.test_vals(gc4_1_low_ru, "timeseries", test=False))
gc4_2=(noramp_fit.test_vals(gc4_2_low_ru, "timeseries", test=False))
gc4_3=(noramp_fit.test_vals(gc4_3_low_ru, "timeseries", test=False))
gc4_1a=(noramp_fit.test_vals(gc4_1_alternate, "timeseries", test=False))
gc4_2a=(noramp_fit.test_vals(gc4_2_alternate, "timeseries", test=False))
gc4_3a=(noramp_fit.test_vals(gc4_3_alternate, "timeseries", test=False))
gc_results=[gc4_1, gc4_2, gc4_3]
gc_alternate=[gc4_1a, gc4_2a, gc4_3a]
nondim_times=noramp_fit.t_nondim(time_results)

gc_results=[np.multiply(gc_results[i], noramp_fit.nd_param.c_I0*1000) for i in range(0,len(gc_results))]
gc_dicta=[np.multiply(gc_alternate[i], noramp_fit.nd_param.c_I0*1000) for i in range(0,len(gc_alternate))]
gc_dict=dict(zip(filenames, gc_results))
gc_alternate=dict(zip(filenames, gc_dicta))
def RMSE(data, sim):
    data=data*1000
    error=1.0/len(data)*np.sqrt(np.sum((np.power(np.subtract(data, sim),2))))
    return error
disped_1=gc4_1
disped_2=gc4_2
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
for i in range(0, len(filenames)):
    #plt.subplot(2,3, i+1)
    plt.plot(voltages[filenames[i]], np.multiply(currents[filenames[i]], 1000),label="Data", color=colors[1])
    plt.plot(voltages[filenames[i]], gc_dict[filenames[i]], label="Low resistance", color=colors[0])
    plt.plot(voltages[filenames[i]], gc_alternate[filenames[i]], label="High Resistance", color=colors[2])
    print filenames[i], RMSE(np.multiply(currents[filenames[i]], 1000), gc_dict[filenames[i]]), RMSE(np.multiply(currents[filenames[i]], 1000), gc_alternate[filenames[i]])
    plt.title(filenames[i])
    plt.xlabel("Voltages(V)")
    plt.ylabel("Current(mA)")
    plt.legend()
    plt.subplot(2,3,i+3+1)
    plt.plot(times[filenames[i]], np.multiply(currents[filenames[i]], 1000), label="Data", color=colors[1])
    plt.plot(times[filenames[i]], gc_dict[filenames[i]], label="Low resistance", color=colors[0])
    plt.plot(times[filenames[i]], gc_alternate[filenames[i]], label="High Resistance", color=colors[2])
    plt.plot(times[filenames[i]], np.subtract(np.multiply(currents[filenames[i]], 1000), gc_dict[filenames[i]]),label="Low Resistance Residual" ,color=colors[3])
    plt.plot(times[filenames[i]], np.subtract(np.multiply(currents[filenames[i]], 1000), gc_alternate[filenames[i]]),label="High Resitance Residual",color=colors[4] )
    plt.ylabel("Current(mA)")
    plt.xlabel("Time(s)")
    plt.legend()
    error=np.sum(np.power(np.subtract(np.multiply(currents[filenames[i]], 1000), gc_dict[filenames[i]]),2))
    print error
plt.show()

disped1_harmonics=harm_class.generate_harmonics(time_results, disped_1)
disped2_harmonics=harm_class.generate_harmonics(time_results, disped_2)
noramp_fit.variable_returner()
#harm_class.plot_harmonics(time_results, method="phased", Experimental=data_harmonics, Ramp_free=exp_harmonics, Ramped=ramp_harmonics)
harm_class.harmonics_and_voltages(noramp_fit.t_nondim(time_results), noramp_fit.e_nondim(voltage_results),folder, "phased", \
                            Experimental_harmonics=noramp_fit.i_nondim(data_harmonics), Experimental_time_series=noramp_fit.i_nondim(current_results),\
                            LowzResistance_harmonics=noramp_fit.i_nondim(disped1_harmonics), \
                            LowzResistance_time_series=noramp_fit.i_nondim(disped_1))

fig, ax=plt.subplots(1, 1)
line_elements=[]
"""
l1,=(ax.plot(voltage_results, disped_1, lw=2))
ax.plot(voltage_results, current_results, lw=2, alpha=0.5)



axcolor = 'lightgoldenrodyellow'
slider_ax=[]
slider_ax_element=[]
for i in range(0, len(noramp_fit.optim_list)):
    slider_ax.append(plt.axes([0.25, 0.0+(i*0.02), 0.65, 0.01], facecolor=axcolor))
    slider_ax_element.append(Slider(slider_ax[i], noramp_fit.optim_list[i], param_bounds[noramp_fit.optim_list[i]][0], param_bounds[noramp_fit.optim_list[i]][1],carbon_means_abs[i]) )



def update(val):
    params=np.zeros(len(noramp_fit.optim_list))
    for i in range(0, len(noramp_fit.optim_list)):
        params[i]=slider_ax_element[i].val
    print list(params)
    test=noramp_fit.test_vals(params, likelihood="timeseries", test=False)
    l1.set_ydata(test)
    fig.canvas.draw_idle()

for i in range(0, len(noramp_fit.optim_list)):
    slider_ax_element[i].on_changed(update)
plt.show()
noramp_fit.optim_list=['E_0','k_0', 'Ru',"Cdl","CdlE1", "CdlE2",'gamma', 'omega', "phase", "cap_phase", "alpha"]
ramp_free_fixed_Ru=(noramp_fit.test_vals(carbon_means_6_20_ru, "timeseries", test=False))

test=ramp_free_fixed_Ru
noise_val=0.00
noise_max=max(test)*noise_val
noise=np.random.normal(0,noise_max, len(test))
synthetic_data=np.add(ramp_free_fixed_Ru, noise)
fourier_test1=noramp_fit.kaiser_filter(synthetic_data)
test_data=np.fft.ifft(likelihood_func)
L=len(test)
time_len=range(0, len(time_results))
f_len=np.linspace(0, time_results[-1], len(fourier_test1))
time_plot=np.interp(f_len, time_len, time_results)
hann=np.hanning(L)
"""
dummy_times=np.linspace(0, 1, len(likelihood_func))
#noramp_fit.optim_list=['Ru', 'omega']
nodisp_results=carbon_means_disped_gc43
noramp_fit.def_optim_list(["E_0", "k_0","Ru","Cdl","CdlE1","CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])
noramp_fit.test_vals(nodisp_results, "timeseries")
noramp_fit.def_optim_list(["E0_mean","E0_std", "k_0","Ru","Cdl","CdlE1", "CdlE2",'gamma', 'omega',"phase","cap_phase","alpha"])
#noramp_fit.def_optim_list(["E0_mean", "E0_std", "Ru", "alpha"])
noramp_fit.dim_dict["CdlE1"]=0
noramp_fit.dim_dict["CdlE2"]=0
noramp_fit.dim_dict["CdlE3"]=0
#noramp_fit.def_optim_list(["E0_mean", "E0_std", "Ru", "k_0"])
fourier_arg=likelihood_func
true_data=current_results
#if simulation_options["experimental_fitting"]==False:d
#elif simulation_options["experimental_fitting"]==True:
    #fourier_arg=likelihood_func
    #true_data=current_results
noramp_fit.pass_extra_data(true_data, fourier_arg)
if simulation_options["likelihood"]=="timeseries":
    cmaes_problem=pints.SingleOutputProblem(noramp_fit, time_results, true_data)
elif simulation_options["likelihood"]=="fourier":
    cmaes_problem=pints.SingleOutputProblem(noramp_fit, dummy_times, fourier_arg)

score = pints.SumOfSquaresError(cmaes_problem)#[4.56725844e-01, 4.44532637e-05, 2.98665132e-01, 2.96752050e-01, 3.03459391e-01]#
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(noramp_fit.optim_list))], [np.ones(len(noramp_fit.optim_list))])
noramp_fit.simulation_options["label"]="cmaes"
x0=abs(np.random.rand(noramp_fit.n_parameters()))
print x0
print len(x0), noramp_fit.n_parameters()
num_runs=1
score_mat=np.zeros(num_runs)
noramp_fit.variable_returner()
param_mat=np.zeros((num_runs,len(noramp_fit.optim_list)))
score_vec=np.zeros(num_runs)
print file
for i in range(0, num_runs):
    x0=abs(np.random.rand(noramp_fit.n_parameters()))#noramp_fit.change_norm_group(gc4_3_low_ru, "norm")
    print noramp_fit.change_norm_group(x0, "un_norm")
    cmaes_fitting=pints.Optimisation(score, x0, sigma0=None, boundaries=CMAES_boundaries, method=pints.CMAES)
    cmaes_fitting.set_max_unchanged_iterations(iterations=200, threshold=1e-3)
    found_parameters, found_value=cmaes_fitting.run()
    cmaes_results=noramp_fit.change_norm_group(found_parameters, "un_norm")
    print list(cmaes_results)
    cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
    #plt.subplot(3,5,i+1)
    print list(cmaes_results)
    #noramp_fit.simulate(found_parameters,time_results, normalise=True, likelihood="fourier", test=True )
    cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
    cmaes_fourier=noramp_fit.test_vals(cmaes_results, likelihood="fourier", test=False)
    param_mat[i,:]=cmaes_results
    score_vec[i]=found_value
    #plt.title(file)
    #plt.subplot(1,2,1)
    plt.plot(voltage_results, true_data)
    plt.plot(voltage_results, cmaes_time)
    #plt.subplot(1,2,2)
    #plt.plot(time_results, true_data)
    #plt.plot(time_results, cmaes_time)
    plt.show()
    #fourier_data=np.fft.ifft(fourier_arg)
    #results=np.fft.ifft(cmaes_fourier)
    #plt.plot(fourier_data)
    #plt.plot(results)
    #plt.show()
print file
print "low_Ru"
for i in range(0, len(param_mat)):
    print i
    print param_mat[i,:]
    print score_vec[i]
    print "---" *10
#best_idx=np.where(score_vec==min(score_vec))
#best_idx=best_idx[0][0]
cmaes_results=gc4_1_alternate#param_mat[i,:]
#cmaes_results=np.array([0.2507124585192858, 0.012810161448688611, 81.25609725526283,100, 3.4455963089647964e-05, 0.08570032724521258, -0.0026762829980299324, 1.4591816196158548e-10, 8.940730783653084, 4.965641139765749, 4.375333040247454, 0.7811192097028208])#[0.2546334543150689, 0.03323126819215362, 116.24098132960108,144.78687238133318,3.4422425270083556e-05, 0.09197677289665342, -0.002971108200729257, 1.4130483196036997e-10,8.94090414236465, 4.947246577563367, 4.361083156674927, 0.5])#param_mat[best_idx, :]
#high alpha disped [0.2507124585192858, 0.012810161448688611, 81.25609725526283,100, 3.4455963089647964e-05, 0.08570032724521258, -0.0026762829980299324, 1.4591816196158548e-10, 8.940730783653084, 4.965641139765749, 4.375333040247454, 0.7811192097028208]
print list(cmaes_results)
cmaes_time=noramp_fit.test_vals(cmaes_results, likelihood="timeseries", test=False)
plt.plot(voltage_results, current_results)
plt.plot(voltage_results, cmaes_time)
plt.show()
#error=np.std(np.subtract(cmaes_prediction, likelihood_func))

error=np.std(np.subtract(cmaes_time, current_results))
mcmc_problem=pints.SingleOutputProblem(noramp_fit, time_results, current_results)
#mcmc_problem=pints.SingleOutputProblem(noramp_fit, time_results, current_results)
#updated_lb=np.append(cmaes_results*0.75, [0])#found_parameters[3]*0.97,
#updated_ub=np.append(cmaes_results*1.25, [2*error])#found_parameters[3]*1.03,
#updated_boundaries=[updated_lb, updated_ub]
#updated_boundaries=np.sort(updated_boundaries, 0)

#noramp_fit.define_boundaries(updated_boundaries)
updated_lb=np.append(noramp_fit.boundaries[0],0)
updated_ub=np.append(noramp_fit.boundaries[1], 2*error)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_lb,
                                updated_ub)
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(cmaes_results, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
noramp_fit.simulation_options["label"]="MCMC"
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()
f=open("GC4_MCMC_1_high_ru", "w")
np.save(f, chains)
f.close()
pints.plot.trace(chains)
plt.show()

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
