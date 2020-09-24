import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm
for j in [0.1,0.25, 0.35, 0.5]:
    size=1000
    k0_weights=np.zeros(size)
    shape=j
    scale=1

    loc=0
    k_start=lognorm.ppf(1e-9, shape, loc, scale=scale)
    k_end=lognorm.ppf(1-1e-9, shape, loc, scale=scale)
    k0_vals=np.linspace(k_start,k_end, size)
    k0_midpoints=np.zeros(len(k0_vals))
    k0_weights[0]=lognorm.cdf(k0_vals[0], shape, loc, scale=scale)
    for k in range(1, len(k0_weights)):
        k0_weights[k]=lognorm.cdf(k0_vals[k], shape, loc, scale=scale)-lognorm.cdf(k0_vals[k-1], shape, loc, scale=scale)
        k0_midpoints[k]=(k0_vals[k]+k0_vals[k-1])/2
    print(k0_midpoints)
    #plt.axvline(j)
    #plt.scatter(np.log(k0_midpoints), np.cumsum(k0_weights), label=j)
    plt.plot(np.log(k0_midpoints),np.cumsum(k0_weights), label=j)
print(sum(k0_weights))
plt.xlabel("log(X)")
plt.ylabel("CDF(X)")
plt.legend()
plt.show()
