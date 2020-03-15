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
coverage=5e-7
e0_val_range=np.linspace(norm.ppf(coverage, loc=location, scale=scale),norm.ppf(1-coverage, loc=location, scale=scale), length)
e0_val_range[0]=norm.ppf(coverage/2, loc=location, scale=scale)
e0[0]=norm.cdf(e0_val_range[0], loc=location, scale=scale)
k0=np.zeros(length)
#plt.plot(x, lognorm.cdf(x, s, loc=0, scale=1),'r-', lw=5, alpha=0.6, label='lognorm pdf')
#plt.show()
start=time.time()
e0_val_range2=np.zeros(len(e0_val_range))

width=np.mean(np.diff(e0_val_range))
midpoints=e0_val_range[:-1] + np.diff(e0_val_range)/2
for i in range(1, len(e0)):
#    e0_val_range2[i]=(e0_val_range[i]+e0_val_range[i-1])/2
    e0[i]=norm.cdf(e0_val_range[i], loc=location, scale=scale)-norm.cdf(e0_val_range[i-1], loc=location, scale=scale)
    #e0[i]=(e0_val_range[i]-e0_val_range[i-1])*norm.pdf(e0_val_range2[i], loc=location, scale=scale)

#plt.plot(e0_val_range2, e0)#, width=0.05)
pdf=norm.pdf(e0_val_range, loc=location, scale=scale)
#pdf=e0
print(sum(e0))
plt.plot(e0_val_range, pdf, color="blue", lw=2)
bins=16
e0_bins=np.linspace(norm.ppf(coverage, loc=location, scale=scale),norm.ppf(1-coverage, loc=location, scale=scale), bins)
pdf_vals=np.interp(e0_bins, e0_val_range, e0)
pdf_vals=norm.pdf(e0_bins, loc=location, scale=scale)

_, ymax = plt.ylim()
line_len=100
current_line=[e0_bins[0]]*line_len
chosen_bin=9
for i in range(0,bins):
    current_line_x=[e0_bins[i]]*line_len
    current_line_y=np.linspace(0, pdf_vals[i], line_len)
    plt.plot(current_line_x,current_line_y, color="black")
    if i==chosen_bin:
        prev_line_x=[e0_bins[i-1]]*line_len
        pdf_idx=tuple(np.where((e0_val_range>prev_line_x[0])& (e0_val_range<current_line_x[0])))
        pdf_y_coords=pdf[pdf_idx]
        pdf_x_coords=e0_val_range[pdf_idx]
        bottom_line_x=pdf_x_coords
        plt.fill_between(bottom_line_x, pdf_y_coords)

#plt.axvline(location)
plt.ylim(ymin=0)
plt.xlabel("X")
plt.ylabel("Probability density")
fig=plt.gcf()
fig.set_size_inches((7, 4.5))
plt.show()
save_path="Dispersion example"
fig.savefig(save_path, dpi=500)
