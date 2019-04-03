import  matplotlib.pyplot as plt
import pints.plot
import numpy as np

length_list=[1e4, 2e4]
dec_list=[8]
repeat_num=20
for lcv_1 in range(0, len(length_list)):
    for lcv_2 in range(0, len(dec_list)):
        for lcv_3 in range(0, repeat_num):
            desired_length=int(length_list[lcv_1])
            dec_amount=dec_list[lcv_2]
            filename=str(desired_length)+"_"+str(dec_amount)+"_"+str(lcv_3)+".Black"
            print filename
            chains=np.load(filename)
            mean_params=np.zeros(len(chains[0, 0,:]))
            for i in range(0, len(chains[0, 0,:])):
                mean_params[i]= np.mean(chains[:,5000:,i])
            print mean_params

            print filename
            pints.plot.trace(chains)
            plt.show()
"""
[2.32335413e-01 5.68179665e+00 4.89371647e+02 9.99490798e-05 2.07607740e-10 1.31008993e+03]
[2.86962985e-01 4.94683058e+00 6.91546540e+02 4.67094486e-05 1.88261079e-10 1.06230370e+03]
[9.22670363e-02 3.20189625e+01 2.98799320e+02 4.63890146e-04 6.02505520e-10 2.79938713e+03]
[2.32347974e-01 5.62703447e+00 4.89017341e+02 9.93453389e-05 2.07207048e-10 2.33081802e+03]
[2.87868599e-01 5.17738747e+00 6.82499495e+02 4.97483777e-05 1.92829232e-10 2.68775418e+03]
[4.49168483e-01 7.93188672e+01 8.12297124e+02 11e-4.11441988e-04 5.41337106e-10 1.18941002e+03]
"""
