import pbwrap.quantification as quant
import numpy as np
import pbwrap.utils as utils
import pbwrap.quantification.cell as cellquant
import numpy as np
import pbwrap.quantification.measures as measures
from skimage.measure import regionprops_table
from scipy.signal import fftconvolve
from pbwrap.quantification.measures import count_spots_in_mask, count_rna_close_pbody_global
import time
import pandas as pd

"""
Testing : 26/10/23 ; ###Attempt at counting number of spots within given radius of a pixel, for all pixels
"""
rng = np.random.default_rng(seed=100)
anchor = np.zeros((30,2048, 2048))
spots = np.zeros((30,2048, 2048))
Z, Y, X = rng.integers(0,10,(3,10000))
anchor_list = list(zip(Z,Y,X))
anchor[Z, Y,X] += 1

Z, Y, X = rng.integers(0,10,(3,10000))
spots_list = list(zip(Z,Y,X))
spots[Z, Y,X] += 1

clock = time.process_time()
count = measures.spots_colocalisation_à_rename(spots_list, anchor_list, radius_nm= 310, image_shape= (30, 2048, 2048), voxel_size= (300,103,103))
print('time : ', time.process_time() - clock)
print(count)


# """

# #spot counting
# mask = np.zeros([4,4])
# mask[1:,0:2] = 1
# print(mask)

# spots = [[0,0], [1,1], [1,2], [2,1], [3,2]]

# res = quant.count_spots_in_mask(spots, mask)
# print(res)
# """

# #signal metrics

# l1 = [0,2,3,4,0]
# l2 = [1,100,1,100,100]
# l3 = [1,1,50,2,1]
# l4 = [1,1,25,3,1]
# l5 = [1,1,2,2,4]

# m1 = [0,0,0,0,0,0,0,0]
# m2 = [0,1,0,1,0,3,3,0]
# m3 = [0,1,1,1,0,0,0,0]
# m4 = [0,0,0,1,0,0,4,0]
# m5 = [2,0,0,0,0,0,4,0]

# data = np.array([l1,l2,l3,l4,l5])
# mask = np.array([m1,m2,m3,m4,m5], dtype= int)

# print('data\n',data)
# print("mask\n", ~mask)

# X = [1,3,1,2,3,3,0, 0,2,6,5]
# Y = [1,1,2,2,2,3,4,0,1,3,3]
# coords = list(zip(Y,X))
# print("coords in mask : \n",mask[Y,X])

# LEFT = pd.DataFrame({'id' : list(np.arange(9)) + [1], 'name' : ['a','b','c','d','e','f','g','h','u','j']})
# RIGHT = pd.DataFrame({'id' : np.random.randint(1,10,10), 'price' : np.random.randint(1,5,10)})

# print(LEFT)
# print(RIGHT)
# print(pd.merge(LEFT,RIGHT, how='left', on= 'id', validate= 'one_to_many'))





#Testing measures.count_rna_close_pbody_global
# distance = [0,1000]#,100,200,400,600,800,1000,1500,2000]
# Pbody_dictionary = count_rna_close_pbody_global(mask, spots_coords= coords, distance_nm= distance, voxel_size=(100,100))
# print(Pbody_dictionary)



# DF = pd.DataFrame({'label' : [1,2,3,4]})
# print(DF)
# for dist in distance :
#     print("Distance : ", dist)
#     print("Res frame : ", Pbody_dictionary['spot {0} nm'.format(dist)])
#     print("Res frame : \n", pd.DataFrame(columns= ['label', 'rna {0}nm count'.format(dist)], data= zip(*Pbody_dictionary['spot {0} nm'.format(dist)])))
#     DF = pd.merge(DF, pd.DataFrame(columns= ['label', 'rna {0}nm count'.format(dist)], data= zip(*Pbody_dictionary['spot {0} nm'.format(dist)])), how= 'left', on= 'label')
#     # DF = pd.merge(DF, pd.DataFrame(columns= ['label', 'malat1 {0}nm count'.format(dist)], data= Pbody_dictionary['malat1 {0} nm']), how= 'left', on= 'label')
# DF = DF.fillna(0)
# print(DF)