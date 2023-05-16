import pbwrap.quantification as quant
import numpy as np
import pbwrap.utils as utils
from pbwrap.quantification.cell import compute_signalmetrics
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
l2 = [1,1,1,1,100]
l3 = [1,1,50,2,1]
l4 = [1,1,25,3,1]
l5 = [1,1,2,2,4]

m1 = [0,0,0,0,0]
m2 = [0,1,0,1,0]
m3 = [0,1,1,1,0]
m4 = [0,0,0,1,0]
m5 = [1,0,0,0,0]

data = np.array([l1,l2,l3,l4,l5])
mask = np.array([m1,m2,m3,m4,m5])


print('data\n',data)
print("mask\n", mask)

truth = utils.point_is_in_mask([2,1],mask)

print(truth)