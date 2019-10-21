import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm
for j in [20**x for x in range(1, 4)]:
    size=j
    k0_weights=np.zeros(size)
    shape=0.5
    scale=10000
    
    loc=0
    k_start=lognorm.ppf(0.0001, shape, loc, scale=scale)
    k_end=lognorm.ppf(0.9999, shape, loc, scale=scale)
    k0_vals=np.linspace(k_start,k_end, size)
    k0_weights[0]=lognorm.cdf(k0_vals[0], shape, loc, scale=scale)
    for k in range(1, len(k0_weights)):
        k0_weights[k]=lognorm.cdf(k0_vals[k], shape, loc, scale=scale)-lognorm.cdf(k0_vals[k-1], shape, loc, scale=scale)

    plt.plot(k0_vals, k0_weights)
    print(sum(k0_weights))
plt.show()
