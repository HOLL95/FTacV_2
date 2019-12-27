from multiplotter import multiplot
import matplotlib.pyplot as plt
plot=multiplot(10, 4, **{"harmonic_position":1, "num_harmonics":5, "orientation":"portrait", "fourier_position":0})
print(plot.axes_dict.keys())
