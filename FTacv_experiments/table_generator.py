import os
from decimal import Decimal
import pickle
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
    'gamma': 'mol\:cm^{-2}$',
    'k_0': 's^{-1}$', #(reaction rate s-1)
    'alpha': "",
    "E0_mean":"V",
    "E0_std": "V",
    "k0_shape":"",
    "k0_loc":"",
    "k0_scale":"",
    "cap_phase":"rads",
    'phase' : "rads",
    "":"",
    "error": "mA"
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
    'phase' : "Phase",
    "":"Experiment",
    "error":"Error",
}
optim_list=["",'E0_mean', "E0_std","k_0","Ru", "Cdl","CdlE1","CdlE2", "gamma","phase", "cap_phase", "alpha"]
name_list=[fancy_names[x] for x in optim_list]
values=[[0.23238700571881385, 0.049551518548422144, 115.05221454733167, 65, 7.618955326361026e-05, 0.008505818667880971, -0.0007155635298508722, 8.768552294174988e-11, 4.220784514175496, 4.8275937272117835, 0.5], [0.20577567177631131, 0.062479806282236386, 128.8459667852348, 195.33266794539094, 1.11443470348505e-05, -0.04999999906315923, 0.009999999940122661, 6.763627615826423e-11, 6.283185307170834, 6.283185306767232, 0.5999999982019539], [0.22571662943838539, 0.02493906225072903, 79.58827176810736, 639.9132040714545, 7.571223157005296e-05, 0.006402229706228689, -0.0005824039425956027, 8.140433631811794e-11, 4.381472254429238, 5.098795186930785, 0.5998761440565236]]
names=["Pinned parameters", "Ramped inference", "Fourier objective"]

#names=my_list[0]
f =open("image_tex_test.tex", "r")
table_file=open("image_tex_edited.tex", "w")
param_num=len(name_list)+1
table_title="\\multicolumn{"+str(param_num-1)+"}{|c|}{Parameter values}\\\\ \n"
table_control1="\\begin{tabular}{|"+(("c|"*param_num))+"}\n"
table_control2="\\begin{tabular}{|"+(("p{3cm}|"*param_num))+"}\n"
value_rows=[]
row_1=""
for i in range(0, len(name_list)):
    if unit_dict[optim_list[i]]!="":
        unit_dict[optim_list[i]]=" ("+unit_dict[optim_list[i]]+")"
    if i ==len(values[0]):
        row_1=row_1+name_list[i]+unit_dict[optim_list[i]] +"\\\\\n"
    else:
        row_1=row_1+(name_list[i]+unit_dict[optim_list[i]]+" & ")
value_rows.append(row_1)
for i in range(0, len(names)):
    row_n=""
    row_n=row_n+(names[i]+ " & ")
    for j in range(0, len(values[0])):
        if j ==len(values[0])-1:
            end_str="\\\\\n"
        else:
            end_str=" & "
        print names[i]
        if abs(values[i][j])>1e-2:

            row_n=row_n+(str(round(values[i][j],3))+ end_str)
        else:
            row_n=row_n+("{:.3E}".format(Decimal(str(values[i][j])))+ end_str)
    value_rows.append(row_n)

control_num=0

for line in f:
    if line[0]=="\\":
        if line.strip()=="\hline":
            table_file.write(line)
        try:
            line.index("{")
            command=line[line.index("{")+1:line.index("}")]
            if (command=="tabular")and (control_num==0):
                line=table_control1
                control_num+=1
            elif (command=="tabular") and (control_num==2):
                line=table_control2
                control_num+=1
            elif command=="4":
                line=table_title
            table_file.write(line)
        except:
            continue
    elif line[0]=="e":
        for q in range(0, len(names)+1):
            line=value_rows[q]
            table_file.write(line)
            table_file.write("\hline\n")
f.close()
table_file.close()
filename=""
filename="alice_yellow_3_params_"
filename=filename+"table.png"
os.system("pdflatex image_tex_edited.tex")
os.system("convert -density 300 -trim image_tex_edited.pdf -quality 100 " + filename)
os.system("mv " +filename+" ~/Documents/Oxford/FTacV_2/FTacV_results/Param_tables")
