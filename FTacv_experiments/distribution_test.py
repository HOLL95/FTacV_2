from scipy.stats import norm, lognorm
import matplotlib.pyplot as plt
import numpy as np
import time
location=10000
scale=100
s=0.1
k_scale=np.exp(1)
k_loc=0.954
length=1000
e0=np.zeros(length)

e0_val_range=np.linspace(7000, 12000, length)
e0[0]=norm.cdf(e0_val_range[0], loc=location, scale=scale)
k0_val_range=np.linspace(0,1e4, length)
k0=np.zeros(length)
x = k0_val_range
#plt.plot(x, lognorm.cdf(x, s, loc=0, scale=1),'r-', lw=5, alpha=0.6, label='lognorm pdf')
#plt.show()
start=time.time()

for i in range(1, len(e0)):
    e0[i]=norm.cdf(e0_val_range[i], loc=location, scale=scale)-norm.cdf(e0_val_range[i-1], loc=location, scale=scale)
plt.plot(e0_val_range, e0)
plt.show()
s=1.0
for i in np.arange(1, 1e4, 1000):
    for j in np.arange(100, 2e3, 100):
        k0[0]=lognorm.cdf(k0_val_range[0], s, loc=j, scale=i)
        for k in range(1, len(k0)):
            k0[k]=lognorm.cdf(k0_val_range[k], s, loc=i, scale=j)#-lognorm.cdf(k0_val_range[k-1], s, loc=i, scale=j)
        print("S", s,"location",i, "scale",j, np.sum(k0))
        plt.plot(k0_val_range, k0, label=np.sum(k0))
    #plt.legend()
    plt.show()
