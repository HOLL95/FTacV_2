from scipy.stats import norm, lognorm
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.hermite import hermgauss
import math
import time

location=10
scale=0.5
s=0.1
length=100
nodes, weights=hermgauss(5)
e0=np.zeros(length)
def generic_func(x):
    return x**3+(2*x+1)
def gh_convert(weight, node, loc, scale):
    return (1/math.sqrt(math.pi))*weight*generic_func((scale*math.sqrt(2)*node)+loc)
e0_val_range=np.linspace(norm.ppf(1e-5, loc=location, scale=scale),norm.ppf(1-1e-5, loc=location, scale=scale), length)
e0[0]=norm.cdf(e0_val_range[0], loc=location, scale=scale)
k0=np.zeros(length)
#plt.plot(x, lognorm.cdf(x, s, loc=0, scale=1),'r-', lw=5, alpha=0.6, label='lognorm pdf')
#plt.show()
start=time.time()
e0_val_range_2=np.zeros(len(e0_val_range))
e0_val_range_2[0]=e0_val_range[0]
for i in range(1, len(e0)):
    e0[i]=norm.cdf(e0_val_range[i], loc=location, scale=scale)-norm.cdf(e0_val_range[i-1], loc=location, scale=scale)
    e0_val_range_2[i]=(e0_val_range[i]+e0_val_range[i-1])/2
print(sum(e0))
nc_disp=np.multiply(e0, [generic_func(x) for x in e0_val_range_2])
plt.plot(e0_val_range_2,  nc_disp)#, width=0.05)
gh_disp=[gh_convert(w_x, n_x, location, scale) for w_x, n_x in zip(weights, nodes)]
converted_nodes=[(scale*math.sqrt(2)*node)+location for node in nodes]
print(np.sum(gh_disp), np.sum(nc_disp))
plt.semilogy(converted_nodes, gh_disp)
plt.show()
