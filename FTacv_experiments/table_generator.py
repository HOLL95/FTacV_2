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
    "error": "$\mu$A"
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
optim_list=["","E0_mean", "E0_std","k_0","Ru","Cdl","CdlE1", "CdlE2","gamma","omega","cap_phase", "phase", "alpha", "error"]
name_list=[fancy_names[x] for x in optim_list]
values=[[0.23794865088573297, 0.012233237509352422, 127.31889407388738, 866.5594838874783, 7.735179022904297e-05, 0.002100658721047255, -0.0003313772993110888, 7.92266323442856e-11, 8.940635097832619, 4.4488399948655255, 5.120897513324314, 0.5999996893552537, 1.46194756457142] ,
[0.23363378517047495, 0.03481010462713212, 125.2418680556816, 630.0772611187369, 7.694171200331618e-05, 0.003209611999861764, -0.0004263185805571494, 7.476933579891946e-11, 8.940473342512918, 4.384420736979492, 5.02832295276801, 0.5999996422725197, 0.9523616806936814] ,
[0.23422376456952138, 0.03260990314447788, 127.84406477439944, 663.1731267606777, 7.6536077935438e-05, 0.0028147360128534457, -0.00040069155469661145, 7.510752378483546e-11, 8.940491556539339, 4.394518618676131, 5.040933894042927, 0.5999999355676289, 0.9887550185838346] ,
[0.23388480903605996, 0.030956721349239304, 125.32474893610775, 682.3612310981944, 7.633158635035894e-05, 0.0025282530103625106, -0.00038943293495081674, 7.525190180668605e-11, 8.940505020536172, 4.399328632878901, 5.05002670496206, 0.5952709809282939, 0.9948215522414872] ,
[0.23251409481227692, 0.03738294965025301, 131.6521158525072, 607.0775808516498, 7.621231601962323e-05, 0.0026921147914714116, -0.00041118732365409173, 7.377445613940177e-11, 8.940470763542377, 4.3777618299742835, 5.012665587646454, 0.5999999998930867, 0.8902094611790935] ,
[0.22717146209373792, 0.04358172945485626, 133.36976450230978, 547.1596593800596, 7.590152453613036e-05, 0.0035752502405714953, -0.0004639004506735312, 7.236875390583806e-11, 8.940464616230182, 4.349188487161576, 4.983465808953585, 0.5999999983569773, 0.9885229039950675] ,
[0.23109347783067546, 0.040826523925818856, 135.4589250447268, 560.3943687177551, 7.612035213575633e-05, 0.002603243569825754, -0.00042124723317499453, 7.267099821808608e-11, 8.940431787150844, 4.365541963219826, 4.991452901786719, 0.5999999993648816, 0.8287397071167669] ,
[0.23088787720198278, 0.0414757167565475, 135.63607798101629, 552.7801493376082, 7.589058106581087e-05, 0.0026746235362428497, -0.0004254788676831059, 7.216845469620117e-11, 8.940458689513898, 4.3621420865641705, 4.985965784062582, 0.599999989782048, 0.8226518586563483] ,
[0.22569328325036095, 0.045062403352648175, 137.31230641230678, 532.2778140615204, 7.551081070710953e-05, 0.0032200614250021548, -0.0004557696937827378, 7.149183693364643e-11, 8.940446483671328, 4.343666918770795, 4.971839471475184, 0.5999997002744768, 0.9655414681595011] ,
[0.2311133473596265, 0.04121633633689883, 135.00011298017264, 563.35080075368, 7.543637810840417e-05, 0.0025134968399891378, -0.0004133051212718425, 7.173747475738812e-11, 8.940476458104767, 4.363628045570221, 4.989057338834849, 0.599999996269962, 0.83121546990775104] ,
]


names=["Scan\ "+str(i) for i in range(2, 11)]

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
filename="Alice_1_9_"
filename=filename+"table.png"
os.system("pdflatex image_tex_edited.tex")
os.system("convert -density 300 -trim image_tex_edited.pdf -quality 100 " + filename)
os.system("mv " +filename+" ~/Documents/Oxford/FTacV_2/FTacV_results/Param_tables")
