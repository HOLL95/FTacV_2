import itertools
import time
import numpy as np
list1=list(range(0,2))
list2=list(range(2, 4))
list3=list(range(4,6))
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
for element3 in list3:
    for element2 in list2:
        for element1 in list1:
            total_list.append([element1, element2, element3])
print(time.time()-start)
start=time.time()
bigg_list=list1+list2+list3
combine_3d=itertools.product(bigg_list, repeat=3)
#print([x for x in combine_3d])
print(time.time()-start)
cheat_global=0
def loop_rec(y, n, new_list, tracking_var):
    #print("n1",n)
    tracking_var+=1
    print(n)
    if n >= 1:
        for x in range(y):
            tracking_var, _=loop_rec(y, n - 1, new_list, tracking_var)
    else:
        new_list.append(tracking_var)

    return tracking_var, new_list
       #print("n3", n)
start=time.time()
print(loop_rec(5, 3, [], 0))
print(time.time()-start)
