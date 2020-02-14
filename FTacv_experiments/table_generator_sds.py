import os
from decimal import Decimal
import pickle
import numpy as np
unit_dict={
    "E_0": "V",
    'E_start': "V", #(starting dc voltage - V)
    'E_reverse': "V",
    'omega':"Hz",#8.88480830076,  #    (frequency Hz)
    'd_E': "V",   #(ac voltage amplitude - V) freq_range[j],#
    'v': '$s^{-1}$',   #       (scan rate s^-1)
    'area': '$cm^{2}$', #(electrode surface area cm^2)
    'Ru': "$\\Omega$",  #     (uncompensated resistance ohms)
    'Cdl': "F", #(capacitance parameters)
    'CdlE1': "",#0.000653657774506,
    'CdlE2': "",#0.000245772700637,
    'CdlE3': "",#1.10053945995e-06,
    'gamma': 'mol cm^{-2}$',
    'k_0': 's^{-1}$', #(reaction rate s-1)
    'alpha': "",
    "E0_mean":"V",
    "E0_std": "V",
    "k0_shape":"",
    "k0_loc":"",
    "k0_scale":"",
    "cap_phase":"rads",
    'phase' : "rads",
    "alpha_mean": "",
    "alpha_std": "",
    "sigma":"",
    "":"",
    "noise":"",
    "error":"$\\mu A$",
}
fancy_names={
    "E_0": '$E^0$',
    'E_start': '$E_{start}$', #(starting dc voltage - V)
    'E_reverse': '$E_{reverse}$',
    'omega':'$\\omega$',#8.88480830076,  #    (frequency Hz)
    'd_E': "$\\Delta E$",   #(ac voltage amplitude - V) freq_range[j],#
    'v': "v",   #       (scan rate s^-1)
    'area': "Area", #(electrode surface area cm^2)
    'Ru': "Ru",  #     (uncompensated resistance ohms)
    'Cdl': "$C_{dl}$", #(capacitance parameters)
    'CdlE1': "$C_{dlE1}$",#0.000653657774506,
    'CdlE2': "$C_{dlE2}$",#0.000245772700637,
    'CdlE3': "$C_{dlE3}$",#1.10053945995e-06,
    'gamma': '$\\Gamma',
    'k_0': '$k_0', #(reaction rate s-1)
    'alpha': "$\\alpha$",
    "E0_mean":"$E^0 \\mu$",
    "E0_std": "$E^0 \\sigma$",
    "cap_phase":"Capacitance phase",
    "k0_shape":"$k^0$ shape",
    "k0_scale":"$k^0$ scale",
    "alpha_mean": "$\\alpha\\mu$",
    "alpha_std": "$\\alpha\\sigma$",
    'phase' : "Phase",
    "sigma":"$\\sigma$",
    "":"Experiment",
    "noise":"$\sigma$",
    "error":"RMSE",
}

optim_list=["","E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega", "cap_phase","phase",  "alpha_std", "sigma"]
name_list=[fancy_names[x] for x in optim_list]
values=[[0.24153910732168946, 0.044603891113748496, 140.47824310452017, 502.33792593039504, 7.647385675474214e-05, 0.0017712378385968135, -0.00038532095477947643, 7.294913544203554e-11, 8.940503412695545, 4.3487637277733215, 4.957850585551708, 0.17757928994786973, 0.12726029731730748],
[0.24121750374766188, 0.044614142937894535, 141.9583866287712, 501.5623795195233, 7.640462015541115e-05, 0.001573448821389732, -0.00037962537521388234, 7.269210703911166e-11, 8.940489587237728, 4.349140090733003, 4.957190783420997, 0.1783266391512056, 0.12550868476122412],
[0.24097546764255365, 0.04481883685433181, 138.24756157958754, 497.36547270896585, 7.616304688764342e-05, 0.0016539056339240717, -0.00038428564081840794, 7.231888565513295e-11, 8.940504264967933, 4.346756890589177, 4.955871426388552, 0.177330114686983, 0.12480040325639935],
[0.24049810162497975, 0.04465310233293925, 134.771285037017, 498.27830472943623, 7.606134642330584e-05, 0.0014339297596182342, -0.00037788902935737126, 7.230414492381035e-11, 8.940506811448596, 4.346040609502586, 4.957886372731034, 0.17878627096840077, 0.12258072575037722],
[0.24117167531813616, 0.044130560446392623, 135.64446572670033, 511.7165123899569, 7.57152316894211e-05, 0.0014490868743300348, -0.0003706986852740553, 7.20449700317367e-11, 8.94052256122038, 4.3493190681991125, 4.962286820156788, 0.17813247128394036, 0.12627179692999121],
]
sds=[[0.0001169269792943754, 0.0008305894462864009, 3.498925512817276, 11.756176892467957, 6.147869780604016e-08, 0.00023581717944887856, 1.152051571538602e-05, 3.304251109582605e-13, 2.3829953668351576e-05, 0.0033652974648172244, 0.005933245757118431, 0.00730104551513345, 0.0012067491290075544] ,
[0.00011423531675517217, 0.0008293563114805797, 3.702132264532235, 11.746678209357857, 6.167899181528225e-08, 0.0002413333228332799, 1.178377354651386e-05, 3.3514455663647447e-13, 2.2288711888353114e-05, 0.003330585722229247, 0.005992564695067563, 0.007369641616946625, 0.001237062919892189] ,
[0.0001168396518608515, 0.0008096553804910855, 3.48603732823611, 11.600422482145332, 6.160591524935496e-08, 0.00023454767024456153, 1.139517656602675e-05, 3.2591578422565465e-13, 2.3562777929477305e-05, 0.003279634460831279, 0.005872824190067626, 0.007174445434441141, 0.0011958787472538398] ,
[0.00011317078192679746, 0.0008244718426788368, 3.400899469677774, 11.730921478807657, 5.926114780019426e-08, 0.00023264328010241846, 1.1367330968500481e-05, 3.315119467205256e-13, 2.268304504198771e-05, 0.0033155430932922613, 0.005963311325516575, 0.006593645285944483, 0.00119435000013714] ,
[0.00011587817787598058, 0.0008593300890708407, 3.6964352783051333, 12.19608990289096, 6.099952331761285e-08, 0.00023640075188081493, 1.152025125816077e-05, 3.3850644186196096e-13, 2.332493512773059e-05, 0.003428292234150651, 0.006154262515036874, 0.007187996019639223, 0.0011934423691867212]]
parameter_orientation="column"
param_num=len(name_list)+1
names=["Exp\ " + str(x) for x in range(6, 11)]
title_num=len(names)+1
table_file=open("image_tex_edited.tex", "w")


#names=my_list[0]
if parameter_orientation =="column":
    f =open("image_tex_test_param_col.tex", "r")
    table_control_1="\\begin{tabular}{|"+(("c|"*title_num))+"}\n"
    titles=["& {0}".format(x) for x in names]
    titles="Parameter "+(" ").join(titles)+"\\\\\n"

    row_headings=[name_list[i]+" ("+unit_dict[optim_list[i]]+") " if unit_dict[optim_list[i]]!="" else name_list[i] for i in range(1, len(optim_list))]
    numerical_rows=[]
    for j in range(0, len(values[0])):
        int_row=""
        for q in range(0, len(names)):
            sd_val=sds[q][j]
            if sd_val>1e-2:
                sd_str="("+str(round(sd_val, 3))+") "
            else:
                sd_str="("+"{:.3E}".format(Decimal(sd_val))+") "
            if values[q][j]>1e-2:
                int_row=int_row+"& "+(str(round(values[q][j],3)))+" "+sd_str
            else:
                int_row=int_row+"& "+"{:.3E}".format(Decimal(str(values[q][j])))+" "+sd_str

        numerical_rows.append(int_row+"\\\\\n\hline\n")
    for i in range(0, len(numerical_rows)):
        numerical_rows[i]=row_headings[i]+numerical_rows[i]
    for line in f:
        if line[0]=="\\":
            if "begin{tabular}" in line:
                table_file.write(table_control_1)
            else:
                table_file.write(line)
        elif "Parameter_line" in line:
                table_file.write(titles)
        elif "Data_line" in line:
            for q in range(0, len(numerical_rows)):
                table_file.write(numerical_rows[q])



f.close()

table_file.close()
filename=""
filename="Alice_ramped_options_"
filename=filename+"table.png"
os.system("pdflatex image_tex_edited.tex")
os.system("convert -density 300 -trim image_tex_edited.pdf -quality 100 " + filename)
os.system("mv " +filename+" ~/Documents/Oxford/FTacV_2/FTacV_results/Param_tables")
