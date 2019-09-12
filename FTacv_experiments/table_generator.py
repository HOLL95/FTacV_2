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
optim_list=["",'E0_mean', "E0_std","k_0","Ru", "Cdl","CdlE1","CdlE2", "gamma", "omega","phase", "cap_phase", "alpha", "error"]
name_list=[fancy_names[x] for x in optim_list]
values =[[0.24517090736930636, 0.056720112426872436, 223.12968601072467, 269.4152026256152, 3.467489815250535e-05, 0.08073629430570951, -0.003038015139181096, 1.4366133260790924e-10, 8.940854696401248, 4.9005338607773, 4.378474411191984, 0.6999999999960397, 0.5630571757698516],\
    [0.24707364221482686, 0.05611157233858666, 215.93921919728885, 261.0012935458973, 3.457664799584686e-05, 0.08731475943652235, -0.00315097559270702, 1.4305526845381745e-10, 8.940962054668862, 4.897135123391333, 4.379405055574097, 0.6999999999984596,0.5579129672799583],\
    [0.24581388949144, 0.056437609561567396, 220.84863259145277, 268.249787363191, 3.465024308699503e-05, 0.08279171693559563, -0.003087944274776373, 1.4514161563358487e-10, 8.940879185090578, 4.900935582894143, 4.379408688250694, 0.6999999997489996,0.5538567481598893 ],\
    [0.2390424106953814, 0.0634126339128772, 59.56331759126466, 1.0000000052336186, 3.421333538191425e-05, 0.09142331002140631, -0.0034273912843145503, 1.3806680869872118e-10, 8.940855654906262, 4.881140108985157, 4.332212053411549, 0.5726231710329381, 0.5630537020667534],\
    [0.23953731576330606, 0.06251048824509087, 54.80881131862505, 1.0, 3.401227622853829e-05, 0.09658746281397658, -0.0035025931533504253, 1.4039829708596723e-10, 8.940959682968415, 4.88941265694472, 4.331036700670979, 0.5725681192157199, 0.5579095099728557],\
    [0.23953952896099746, 0.06432508696935929, 60.84131849364619, 1.0, 3.42647900201487e-05, 0.08893011038618533, -0.0033602097013514043, 1.3972124127562154e-10, 8.94084399186188, 4.8779396488202185, 4.331377454498156, 0.5724578682161762,0.5538581095139239]]


names=["GC4\ 1\ High\ Resistance", "GC4\ 2\ High\ Resistance", "GC4\ 3\ High\ Resistance", "GC4\ 1\ Low\ Resistance", "GC4\ 2\ Low\ Resistance", "GC4\ 3\ Low\ Resistance"]
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
