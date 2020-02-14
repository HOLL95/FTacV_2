from scipy.stats import norm, lognorm
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.hermite import hermgauss
import time
print(hermgauss(5))
location=0
scale=1
s=0.1
length=10000
e0=np.zeros(length)

e0_val_range=np.linspace(norm.ppf(1e-5, loc=location, scale=scale),norm.ppf(1-1e-5, loc=location, scale=scale), length)
e0[0]=norm.pdf(e0_val_range[0], loc=location, scale=scale)
k0=np.zeros(length)
#plt.plot(x, lognorm.cdf(x, s, loc=0, scale=1),'r-', lw=5, alpha=0.6, label='lognorm pdf')
#plt.show()
start=time.time()
e0_val_range2=np.zeros(len(e0_val_range))
for i in range(1, len(e0)):
    e0[i]=norm.pdf(e0_val_range[i], loc=location, scale=scale)#-norm.pdf(e0_val_range[i-1], loc=location, scale=scale)
    e0_val_range2[i]=(e0_val_range[i]+e0_val_range[i-1])/2
print(sum(e0))
plt.plot(e0_val_range, e0)#, width=0.05)
plt.axvline(location)
plt.show()
