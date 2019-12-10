from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import matplotlib.pyplot as plt
import numpy as np
problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359]]
}

param_values = saltelli.sample(problem, 1)
Y = np.zeros([param_values.shape[0]])
print(param_values)
Y = Ishigami.evaluate(param_values)
plt.plot(Y)
plt.show()
print(Y)
Si = sobol.analyze(problem, Y, print_to_console=True)
