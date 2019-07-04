from scipy.stats import norm, lognorm
import matplotlib.pyplot as plt
import numpy as np
import time
location=0.1
scale=0.05
s=0.1
k_scale=np.exp(1)
k_loc=0.954
length=16
e0=np.zeros(length)

e0_val_range=np.linspace(-0.4, 0.4, length)
e0[0]=norm.cdf(e0_val_range[0], loc=location, scale=scale)
k0_val_range=np.linspace(0,1e4, length)
k0=np.zeros(length)
x = np.linspace(0,1000,length)
#plt.plot(x, lognorm.cdf(x, s, loc=0, scale=1),'r-', lw=5, alpha=0.6, label='lognorm pdf')
#plt.show()
start=time.time()
k0[0]=lognorm.cdf(k0_val_range[0], s, loc=k_loc, scale=k_scale)
for i in range(1, len(e0)):
    e0[i]=norm.cdf(e0_val_range[i], loc=location, scale=scale)-norm.cdf(e0_val_range[i-1], loc=location, scale=scale)
i=0.1
for k in range(1, len(k0)):
    k0[k]=lognorm.cdf(x[k], i, loc=100, scale=200)-lognorm.cdf(x[k-1], i, loc=100, scale=200)
plt.plot(k0_val_range, k0)
plt.show()
e0_mat, k0_mat=np.meshgrid(e0_val_range, k0_val_range)

weight_matrix=np.multiply(e0_mat, k0_mat)
weight_mat=np.zeros((length, length))
print e0_mat, k0_mat
print time.time()-start
print np.sum(k0)
print np.sum(e0)
bins=length
start= time.time()
print time.time()-start
