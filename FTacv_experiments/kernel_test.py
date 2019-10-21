import numpy as np
import matplotlib.pyplot as plt
data_1, data_2=np.random.multivariate_normal([0,1], [[1,0], [0, 100]], 1000).T
xmin=min(data_1)
xmax=max(data_1)
ymin=min(data_2)
ymax=max(data_2)
xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
print(positions)
