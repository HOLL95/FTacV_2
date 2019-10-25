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
optim_list=["","Ru","Cdl","CdlE1", "CdlE2","CdlE3",'omega', "phase"]
name_list=[fancy_names[x] for x in optim_list]
values=[[4.935278186731851, 6.53922500129266e-05, -0.010027952214194474, -0.000120463926091368, 1.6440371959129674e-05, 8.828862457999545, 5.846494748996116],
            [1735.848562286573, 7.387153313864618e-05, -0.010118667659211301, -0.00027396481928436064, 2.3648787718632755e-05, 8.828834554751786, 0]]
names=["High resistance", "Low resistance"]

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
        print(names[i])
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
filename="alice_blank_3_params_"
filename=filename+"table.png"
os.system("pdflatex image_tex_edited.tex")
os.system("convert -density 300 -trim image_tex_edited.pdf -quality 100 " + filename)
os.system("mv " +filename+" ~/Documents/Oxford/FTacV_2/FTacV_results/Param_tables")
