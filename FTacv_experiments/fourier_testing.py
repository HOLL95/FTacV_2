import numpy as np
import matplotlib.pyplot as plt
import math
x=np.linspace(0, 8*math.pi, 1000)
f=np.fft.fftfreq(len(x), x[1]-x[0])
hann=np.hanning(len(x))
amps=[10.0**i for i in np.arange(-3, 1, 1)]
amps=np.flip(amps)
for i in amps:
    y=i*np.sin(x)
    z=np.exp(y)
    Y=np.fft.fft(z)
    plt.semilogy(x, z)
plt.show()
