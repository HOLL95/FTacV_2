import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm
for j in [0.1]:
    size=16
    k0_weights=np.zeros(size)
    shape=0.495
    scale=440.487

    loc=0
    k_start=lognorm.ppf(1e-9, shape, loc, scale=scale)
    k_end=lognorm.ppf(1-1e-5, shape, loc, scale=scale)
    k0_vals=np.linspace(k_start,k_end, size)
    k0_weights[0]=lognorm.cdf(k0_vals[0], shape, loc, scale=scale)
    for k in range(1, len(k0_weights)):
        k0_weights[k]=lognorm.cdf(k0_vals[k], shape, loc, scale=scale)-lognorm.cdf(k0_vals[k-1], shape, loc, scale=scale)
    #plt.axvline(j)
    plt.plot(k0_vals, k0_weights, label=j)
print(sum(k0_weights))
plt.legend()
plt.show()
