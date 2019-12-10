import sympy as sym
import time
import numpy as np
from sympy.utilities.lambdify import lambdify, implemented_function
from sympy import Function
from sympy.abc import x
import mpmath
import math
sym_params=["K", "r", "t", "p0"]
sym_vars=[sym.Symbol(x) for x in sym_params]
variables=dict(zip(sym_params, sym_vars))
Func=variables["K"]/(1+((variables["K"]/variables["p0"])-1)*sym.exp(-variables["r"]*variables["t"]))
derivatives=[sym.diff(Func, x) for x in sym_params]
[print(x,"\n", y,"\n") for x,y in zip(sym_params, derivatives)]
x=sym.symbols("X")
func=x**2
start=time.time()
for i in range(0, 1000):
    func.evalf(subs={"X":i})
print(time.time()-start)
start=time.time()
for i in range(0, 1000):
    np.power(i, 2)
print(time.time()-start)
print(func)
x, y = sym.symbols('x y')
expr = 3*x**2 + sym.log(x**2 + y**2 + 1)
g = sym.lambdify([x], func, modules=['math'])
start=time.time()
for i in range(0, 1000):
    g(i)
print(time.time()-start)
