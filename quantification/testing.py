"""
Testing : 02/02/2024 : scipy.stats testing
"""
import pandas as pd
import numpy as np
import scipy.stats as scstats
import pbwrap.quantification.statistical_test as stats
import CustomPandasFramework.centrosome.analysis as analysis


path = '/home/flo/Downloads/centrosome_scp/'
measure = ['proportion_rna_centrosome', 'index_median_distance_centrosome']
m1,m2 = measure

data = analysis.prepare_data_frame(path, measure, index_keys=[])
data = data.groupby(['rna','treatment'])[m1].apply(list)

enlisted_data = data.reset_index().groupby(['rna'])[m1].apply(list)
print("DATA\n",data)

for rna, rna_data in zip(enlisted_data.index, enlisted_data) :
    print('\n' + rna)
    print('group : {0}'.format(data[rna].index))
    print('ANOVA')
    f_stat,p_value = stats.ANOVA(rna_data)
    print("F-Statistic : {0} \np-value : {1}\n".format(f_stat, p_value))

    print('Tukey-hsd')
    f_stat,p_value = stats.Tukey_hsd(rna_data)
    print("F-Statistic : {0} \np-value : {1}\n".format(f_stat, p_value))

    print("variance")
    var = [np.var(sample) for sample in rna_data]
    print(var)



"""
Testing : 24/01/24 : Closest distance spot1 to spot 2 array

import testing_samples as samp
import pbwrap.quantification.measures as measures

shape = (15,2000,2000)
voxel_size = (300,103,103)
spot1_number = 6000
spot2_number = 10000

gen = samp.get_random_gen(seed= 1)
spot_1 = samp.random_spots_list(spot1_number, shape, gen=gen)
spot_2 = samp.random_spots_list(spot2_number, shape=shape, gen=gen)

signal_1 = measures._reconstruct_spot_signal(shape, spot_1)
signal_2 = measures._reconstruct_spot_signal(shape, spot_2)

distance = measures.closest_spot_distance(spot_1, spot_2, shape, voxel_size=voxel_size)
print(distance)

"""

# """
# Testing : 26/10/23 ; ###Attempt at counting number of spots within given radius of a pixel, for all pixels
# """
# rng = np.random.default_rng(seed=100)
# anchor = np.zeros((30,2048, 2048))
# spots = np.zeros((30,2048, 2048))
# Z, Y, X = rng.integers(0,10,(3,10000))
# anchor_list = list(zip(Z,Y,X))
# anchor[Z, Y,X] += 1

# Z, Y, X = rng.integers(0,10,(3,10000))
# spots_list = list(zip(Z,Y,X))
# spots[Z, Y,X] += 1

# clock = time.process_time()
# count = measures.spots_colocalisation_à_rename(spots_list, anchor_list, radius_nm= 310, image_shape= (30, 2048, 2048), voxel_size= (300,103,103))
# print('time : ', time.process_time() - clock)
# print(count)


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