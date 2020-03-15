import numpy as np
import matplotlib.pyplot as plt

final=10000
a=np.zeros(final)
log_lambda=0
a[0]=log_lambda
t=range(1, final)
adaptations=0

for i in t:
    random=np.random.rand()
    if random>:
        accept=1
        adaptations+=1
    else:
        accept=0

    #print(a[i-1],np.log(a[i-1]), ((i**(-0.6))*(accept-0.25)))
    log_lambda+=np.power(adaptations+1,-0.6)*(accept-0.25)
    a[i]=np.exp(log_lambda)
plot_t=range(0,final)
plt.plot(plot_t,a)
plt.show()
