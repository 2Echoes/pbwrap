import pbwrap.quantification as quant
import numpy as np
import pbwrap.utils as utils
import pbwrap.quantification.cell as cellquant
import numpy as np
from skimage.measure import regionprops_table
from pbwrap.quantification.measures import count_spots_in_mask

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

print('data\n',data)
print("mask\n", ~mask)

X = [1,3,1,2,3,3,0, 0, 2]
Y = [1,1,2,2,2,3,4,0, 1]
print("coords in mask : \n",mask[Y,X])

coords = list(zip(Y,X))
coords_copy = coords.copy()
print("coords :\n",coords)

for point in coords : 
    if not mask[point[0],point[1]] : coords.remove(point)

print("nbre de coords : ", len(coords))
print("nbre de coords_copy : ", count_spots_in_mask(coords_copy, mask))

X_array = np.array(X, dtype= int)[mask[Y,X]]
Y_array = np.array(Y, dtype= int)[mask[Y,X]]

print(X_array)
print(Y_array)
print(np.array(list(zip(Y_array, X_array))))