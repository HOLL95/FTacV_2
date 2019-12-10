
variables=["Cdl", "dE", "Ru", "dI", "gamma", "k0", "theta", "alpha", "E", "I","E0"]
sensitive_params=["Cdl", "Ru", "gamma", "k0", "alpha", "E0"]
symbolic_vars=[sym.Symbol(x) for x in variables]
symbol_dict=dict(zip(variables, symbolic_vars))
d_theta_1=(1-symbol_dict["theta"])*symbol_dict["k0"]*sym.exp((1-symbol_dict["alpha"])*(symbol_dict["E"]-symbol_dict["Ru"]*symbol_dict["I"]-symbol_dict["E0"]))
d_theta_2=(symbol_dict["theta"])*symbol_dict["k0"]*sym.exp(-symbol_dict["alpha"]*(symbol_dict["E"]-symbol_dict["Ru"]*symbol_dict["I"]-symbol_dict["E0"]))
d_theta=d_theta_1+d_theta_2
cap=symbol_dict["Cdl"]*(symbol_dict["dE"]-symbol_dict["Ru"]*symbol_dict["dI"])
didt=(cap+(symbol_dict["gamma"]*d_theta)-symbol_dict["I"])/(symbol_dict["Cdl"]*symbol_dict["Ru"])
didt=sym.simplify(didt)
for i in range(0, len(sensitive_params)):
    derivative=sym.diff(didt, sensitive_params[i])
    print(sym.simplify(sensitive_params[i]))
    print(derivative, "\n")
