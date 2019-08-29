import os
from decimal import Decimal
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
    'gamma': '\\frac{mol}{cm^{2}}$',
    'k_0': 's^{-1}$', #(reaction rate s-1)
    'alpha': "",
    "E0_mean":"V",
    "E0_std": "V",
    "k0_shape":"",
    "k0_loc":"",
    "k0_scale":"",
    "cap_phase":"rads",
    'phase' : "rads",
    "":""
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
    "":"Experiment"
}
optim_list=["",'E0_mean', "E0_std","k_0","Ru", "Cdl","CdlE1","CdlE2", "gamma", "omega","phase", "cap_phase", "alpha"]
name_list=[fancy_names[x] for x in optim_list]
#values=[[0.2507124585192858, 0.012810161448688611, 81.25609725526283,100.0, 3.4455963089647964e-05, 0.08570032724521258, -0.0026762829980299324, 1.4591816196158548e-10, 8.940730783653084, 4.965641139765749, 4.375333040247454, 0.7811192097028208],\
#        [0.2546334543150689, 0.03323126819215362, 116.24098132960108,144.78687238133318,3.4422425270083556e-05, 0.09197677289665342, -0.002971108200729257, 1.4130483196036997e-10,8.94090414236465, 4.947246577563367, 4.361083156674927, 0.5],\
#        [0.24042307692307694, 0.05995192307692309, 81.25609725526283, 100.0, 3.0169689037255904e-05, 0.020629272483467714, -0.0002036780323021597, 1.4591816196158548e-10, 8.959036763806955, 0.0,0.0, 0.5990384615384616],\
#        [0.24032051282051284, 0.0625721153846154, 116.24098132960108,144.78687238133318, 3.0169689037255904e-05, 0.020629272483467714, -0.0002036780323021597, 1.4130483196036997e-10, 8.959036763806955, 0.0,0.0, 0.5],\
#        ]
values=[[0.2507124585192858, 0.012810161448688611, 81.25609725526283,100, 3.4455963089647964e-05, 0.08570032724521258, -0.0026762829980299324, 1.4591816196158548e-10, 8.940730783653084, 4.965641139765749, 4.375333040247454, 0.7811192097028208],\
        [0.25150277022331746, 0.014115779656223786, 78.20813556543256,58.804257252122355, 3.501508623604994e-05, 0.08815081354697163, -0.0030182380193618527, 1.3910336262835592e-10, 8.941274817456382, 4.946748053531174, 4.3670012650213295, 0.7663320409205243],\
        [0.2508929602137271, 0.014962342550521588,  74.5650170269896,13.194027307777567, 3.51489903928365e-05, 0.0898311438945581, -0.0031909694983879135, 1.3515246782381657e-10, 8.94090414236465, 4.933768379207287, 4.358597908430476, 0.7621639225400778]
        ]

names=["GC4\ 1", "GC4\ 2", "GC4\ 3"]
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
for i in range(0, len(names)):
    filename=filename+names[i]+"_"
filename=filename+"table.png"
os.system("pdflatex image_tex_edited.tex")
os.system("convert -density 300 -trim image_tex_edited.pdf -quality 100 " + filename)
os.system("mv " +filename+" ~/Documents/Oxford/FTacV_2/FTacV_results/Param_tables")
