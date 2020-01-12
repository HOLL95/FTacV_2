import itertools
import time
import numpy as np
list1=list(range(0,2))
list2=list(range(2, 4))
list3=list(range(4,8))
list_of_lists=[list1, list2, list3]
start=time.time()
combine_3d=itertools.combinations(list1+list2+list3, 3)

for current_tuple in combine_3d:
    num_elements=[len(set(current_tuple).intersection(set(list))) for list in list_of_lists]
    if all(x==1 for x in num_elements):
        intersect_idxs=[list.index(x) for x, list in zip(current_tuple, list_of_lists)]
print(time.time()-start)
#a=[[list(x)+[y] for x in combine_2d] for y in list3]
#print(a)
start=time.time()
combine_2d=itertools.combinations(list1+list2,2)
total_list=[]
for element in list3:
    total_list.append([list(tuple)+[element] for tuple in combine_2d])
print(time.time()-start)
print(total_list)
