import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
files=os.listdir('.')
#f=open(filename, "r")
titles=['E_0', 'k_0',  'Ru','Cdl']
true_params=[-0.4, 10, 20.0,0.000134,]
noise_vals=[0.01, 0.02, 0.005, 0]
fig, ax=plt.subplots(2,3)
counter=-1
extension=".cm"
chain_select=1
j=0
for j in range(0, len(files)):
    if files[j][-len(extension):-1]== extension[:-1]:
        print((files[j]))
        counter+=1
        chains=np.load(files[j])
        for i in range(0, len(titles)):
            axes=plt.subplot(2,2,i+1)
            labels=files[j][:files[j].index(extension[1])]
            if titles[i]=="sigma":
                chains=np.divide(chains[chain_select,8000:,i], 200*21.921103406)
                axes.hist(chains, alpha=0.4,bins=20, label=labels[:-1])
            else:
                axes.hist(chains[chain_select,5000:,i], alpha=0.4,bins=20, label=labels[:-1])
            axes.set_title(titles[i])
            if titles[i]!='sigma':
                axes.axvline(true_params[i],color="black", linestyle="--")
            axes.legend()
            #axes.set_xlim((0.9*true_params[i], 1.1*true_params[i]))
#axes=plt.subplot(2,3,i+2)
#axes.axis('off')
plt.show()
