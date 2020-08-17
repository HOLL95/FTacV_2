import scipy.special as bessel
import numpy as np
import matplotlib.pyplot as plt
length=10
sum_results=np.zeros(length)
summation=0
for i in range(0,  length):
    modified_val=bessel.iv(i, 0.5*300e-3)
    print(modified_val)
    summation+=modified_val
    sum_results[i]=summation
plt.plot(sum_results)
plt.show()
