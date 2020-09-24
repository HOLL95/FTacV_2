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
values=[
[0.2413977449390575, 0.048182464697725906, 140.03375099744366, 495.7823625045662, 7.612761438466669e-05, 0.0025341144138490705, -0.00042221556534881657, 7.27134356267899e-11, 8.940506405743712, 4.329115987271387, 4.950006457695003, 0.19159097098443717, 0.7533681241433168] ,
[0.24102944844291002, 0.048497349529395314, 141.82521896400752, 491.19384802208276, 7.601409789962726e-05, 0.002477827799330434, -0.0004228455838537012, 7.241395948917027e-11, 8.940493047439322, 4.328279472666, 4.9473156823531115, 0.1921949891403852, 0.7455277540844046] ,
[0.2408140704702503, 0.048202538218362886, 136.38375219019173, 493.49908184684307, 7.581450988788631e-05, 0.0023923283337832085, -0.00041955420136924863, 7.220019074381552e-11, 8.94050702016123, 4.3277735987254475, 4.950151770137567, 0.19115510367237332, 0.7422461522760112] ,
[0.2404444443022998, 0.04782449086187796, 140.48761290690587, 488.04166270023086, 7.583575045270782e-05, 0.001952651694737962, -0.0004025677013963219, 7.200961922103972e-11, 8.940505305256325, 4.331403849422406, 4.945834754904772, 0.18938421569006927, 0.7270175284650521] ,
[0.24103532681889958, 0.04725612528058502, 138.6414711102102, 499.9842386905104, 7.542592287627759e-05, 0.0023002023852481154, -0.0004113355059251488, 7.152736642108646e-11, 8.940522007228264, 4.334097839269449, 4.952057295870433, 0.17471723544230158, 0.7546814619750335] ,
]
sds=[
[5.559697750485778e-05, 0.00036158130628566597, 1.807290719463265, 5.350972196180232, 2.869187248854773e-08, 0.00010157682889483346, 4.964114972160335e-06, 1.4735579298033465e-13, 1.1982733763938402e-05, 0.0015260101789174117, 0.0027953151113136797, 0.003057934459040198, 0.0036053800653138126] ,
[5.55423278884309e-05, 0.0003439049570126469, 1.8408352675188622, 5.0873335013841, 2.7341831249691224e-08, 0.00010379073656422832, 5.113305215798263e-06, 1.4208171194881497e-13, 1.1125147264904768e-05, 0.0014484283179108802, 0.0026545466159533185, 0.0030449205064008804, 0.0035572286052124496] ,
[5.415890024417677e-05, 0.0003500389380409878, 1.7751893416719369, 5.1854947924271375, 2.768332018666361e-08, 0.00010037961540712096, 4.9342882296127384e-06, 1.4156036656553375e-13, 1.1425877835319958e-05, 0.0014695321240623452, 0.0027053757481481895, 0.002890258694849384, 0.003543115030890644] ,
[7.414519958406201e-05, 0.0004956964219279756, 2.5704808709600955, 7.163277897243989, 3.8999547562778335e-08, 0.00014731530000115748, 7.1927895004385755e-06, 1.986566044394388e-13, 1.539918548478229e-05, 0.002018608377087118, 0.003755136825799835, 0.004379982269689827, 0.005119358837658816] ,
[8.230405708724956e-05, 0.0005876740199346263, 2.6500717620882748, 8.792785663763151, 4.2325741805998355e-08, 0.0001624036383442179, 7.849037047140074e-06, 2.328130675240272e-13, 1.6320851157131392e-05, 0.0024724481958989597, 0.00443423957504069, 0.004722254831667136, 0.005154064845321515] ,
]

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
