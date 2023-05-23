import pbwrap.quantification as quant
import numpy as np
import pbwrap.utils as utils
import pbwrap.quantification.cell as cellquant
import time as t
"""

#spot counting
mask = np.zeros([4,4])
mask[1:,0:2] = 1
print(mask)

spots = [[0,0], [1,1], [1,2], [2,1], [3,2]]

res = quant.count_spots_in_mask(spots, mask)
print(res)
"""

#signal metrics

l1 = [0,2,3,4,0]
l2 = [1,100,1,100,100]
l3 = [1,1,50,2,1]
l4 = [1,1,25,3,1]
l5 = [1,1,2,2,4]

m1 = [0,0,0,0,0]
m2 = [0,1,0,1,0]
m3 = [0,1,1,1,0]
m4 = [0,0,0,1,0]
m5 = [1,0,0,0,0]

data = np.array([l1,l2,l3,l4,l5])
mask = np.array([m1,m2,m3,m4,m5], dtype= bool)
mask = np.random.randint(0,2,[5000,5000], dtype= bool)

print('data\n',data)
print("mask\n", ~mask)

timestamp = t.process_time()
#nouvelle solution
pixel_number = mask.sum()
print("Solution1\npixel number = ",pixel_number,"\nProcess time = ",str(t.process_time() - timestamp))

#solution en cours
timestamp = t.process_time()
_, count = np.unique(mask, return_counts= True)
pixel_number = count[1]

print("Solution2\npixel number = ",pixel_number,"\nProcess time = ",str(t.process_time() - timestamp))