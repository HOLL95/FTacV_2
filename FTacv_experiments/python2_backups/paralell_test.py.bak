import numpy as np
from time import time
import multiprocessing as mp
import itertools

"""

def howmany_within_range2(i, row, minimum, maximum):

    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    print i
    return (i, count)

# Prepare data
np.random.RandomState(100)
arr = np.random.randint(0, 10, size=[200, 5])
data = arr.tolist()
print data
data[:5]
pool = mp.Pool(mp.cpu_count())

results = []
def collect_result(result):
    global results
    results.append(result)
for i, row in enumerate(data):
    print i, row
    pool.apply_async(howmany_within_range2, args=(i, row, 4, 8), callback=collect_result)

# Step 4: Close Pool and let all the processes complete
pool.close()
pool.join()  # postpones the execution of next line of code until all processes in the queue are done.
print results
# Step 5: Sort results [OPTIONAL]
results.sort(key=lambda x: x[0])
print results
results_final = [r for i, r in results]

print(results_final[:10])
#Generate values for each parameter
a = range(10)
b = range(10)
c = range(10)
d = range(10)

#Generate a list of tuples where each tuple is a combination of parameters.
#The list will contain all possible combinations of parameters.
paramlist = list(itertools.product(a,b,c,d))

#A function which will process a tuple of parameters
def func(params):
  a = params[0]
  b = params[1]
  c = params[2]
  d = params[3]
  return a*b*c*d

#Generate processes equal to the number of cores
pool = multiprocessing.Pool()

#Distribute the parameter sets evenly across the cores
res  = pool.map(func,paramlist)
#> [3, 1, 4, 4, 4, 2, 1, 1, 3, 3]
"""
class someClass(object):
   def __init__(self, y):
       self.y=y
   def f(self, x):
       return x*x

   def go(self):
      p = mp.Pool(4)
      for i in range(4):
        sc = p.apply_async(self, args=i)
      print sc

   def __call__(self, x):
     return self.f(x)

sc = someClass(5)
sc.go()
