"""
Testing : 02/02/2024 : scipy.stats testing
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as scstats
import pbwrap.plot.stats as plot
import pbwrap.quantification.statistical_test as stats
import CustomPandasFramework.centrosome.analysis as analysis


path = '/home/flo/Downloads/centrosome_scp/'
measure = ['proportion_rna_centrosome', 'index_median_distance_centrosome']
m1,m2 = measure

analysis.p_value(path, path, measure)

# data = analysis.prepare_data_frame(path, measure, index_keys=[])
# data = data.groupby(['rna','treatment'])[m1].apply(list)

# enlisted_data = data.reset_index().groupby(['rna'])[m1].apply(list)
# print("DATA\n",data)
# print(data.index.levels)

# plt.figure(figsize= (15,15))
# ax = plt.gca()
# ax = plot.ANOVA_Tukey_hsd(ax, data, h0_treatment= 'DMSO')
# plt.show()

# for rna, rna_data in zip(enlisted_data.index, enlisted_data) :
#     print('\n' + rna)
#     print('group : {0}'.format(data[rna].index))
#     print('ANOVA')
#     f_stat,p_value = stats.ANOVA(rna_data, return_f_stats=True)
#     print("F-Statistic : {0} \np-value : {1}\n".format(f_stat, p_value))

#     print('Tukey-hsd')
#     f_stat,p_value = stats.Tukey_hsd(rna_data, return_f_stats=True)
#     print("F-Statistic : {0} \np-value : {1}\n".format(f_stat, p_value))

#     print("variance")
#     var = [np.var(sample) for sample in rna_data]
#     print(var)

# plt.figure(figsize= (15,15))
# ax = plt.subplot(2,1,2)
# plot.variance_plot(ax=ax, group_list= enlisted_data, labels= enlisted_data.index, normalise_var=True)
# ax = plt.subplot(2,2,1)
# plot.variance_plot(ax=ax, group_list= enlisted_data, labels= enlisted_data.index, normalise_var=True)
# ax = plt.subplot(2,2,2)
# plot.variance_plot(ax=ax, group_list= enlisted_data, labels= enlisted_data.index, normalise_var=True)
# plt.show()

"""### testing colocalisation_plot ###
import pbwrap.plot.visuals as vis
import testing_samples as samp
import bigfish.stack as stack

seed= 1
gen = samp.get_random_gen(seed=seed)
shape = (15,1000,1000)
spots1_number = 10000
spots2_number = 25
path = '/home/floricslimani/Documents/testing.tif'

spots1 = samp.random_spots_list(spots1_number, shape = shape, gen=gen)
spots2 = samp.random_spots_list(spots2_number, shape = shape, gen=gen)

signal = samp.random_spot_signal(shape, 100)

vis.colocalisation_plot(signal, shape, (1,1,1), 10, path, spots1, spots2, spots2)
read = stack.read_image(path)
print(read.shape)


"""


############
# in_path = "/home/floricslimani/Documents/Projets/1_P_body/stack_O8_p21/output/20230531 17-01-21/result_tables"
# output_path = "/home/floricslimani/Documents/Projets/1_P_body/Workshop/"
# if not in_path.endswith('/') : in_path += '/'
# if not output_path.endswith('/') : output_path += '/'

# # Cleaning tables
# Acquisition = pd.read_feather(in_path + 'Acquisition')
# Cell = pd.read_feather(in_path + 'Cell')
# rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(Acquisition, Cell, 80)
# rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(rna_clean_Acquisition, rna_clean_Cell, 250, keep='lower')
# thresholds = list(Acquisition.loc[:, "RNA spot threshold"])
# gene_list = list(Acquisition.loc[:, "rna name"])
# malat_clean_Acquisition, malat_clean_Cell = update.from_malat_remove_acquisition(Acquisition, Cell, limit= 10)
# new_Cell = update.from_nucleus_malat_proportion_compute_CellullarCycleGroup(malat_clean_Cell)



# # X = np.concatenate([np.arange(100), np.arange(100,0,-1)])
# # Y = (X*np.random.randint(0,100,200)) **2 + 24
# # slope,intercept = simple_linear_regression(X,Y)
# # plt.plot(X,Y,'k.', label = 'random distrib')
# # plt.plot(X,slope*X+intercept, 'r', label = "{0}x + {1}".format(slope,intercept))
# # plt.axis('tight')
# # plt.legend()
# # plt.show()



# scatter.DapiSignal_vs_CellNumber(Cell, integrated_signal= False)


# # scatter.Malat_inNuc_asDapiIntensity(Acquisition= Acquisition.query("`rna name` in ['NF1']"), Cell= Cell, plot_linear_regression= True, title= 'NF1', s= 12)
# # scatter.Malat_inNuc_asDapiIntensity(Acquisition= Acquisition.query("`rna name` in ['PABPC1']"), Cell= Cell, plot_linear_regression= True, title= 'PABPC1', s=12)



# # box.box_plot(new_Cell.loc[:, "malat1 spots in nucleus"])


# #Testing overlapping annotations

# import matplotlib.pyplot as plt
# import numpy as np

# fig = plt.figure(figsize=(10,10))

# #Data
# X = [0,1,2,3,4,5,5.1,5.15,5.12,5.23,5.4,6,7,7.2,8,9,10]
# Y = X

# plt.scatter(X,Y)
# text_post = list(zip(X,Y))
# annotations = ['GENE {0}'.format(i) for i in range(1,len(text_post)+1)]
# obj_list = []
# # for text,pos in zip(annotations,text_post) :
# #     obj_list.append(plt.annotate(text, pos))

# fig.canvas.draw()
# plt.tight_layout()
# xmin, xmax, ymin, ymax = plt.axis()
# xmin = 0
# ymin = 0
# plt.axis([xmin,xmax,ymin,ymax])

# # plt.show()

# #Testing

# pos_list = list(zip(X,Y))
# text_list = annotations
# x_unit, y_unit = compute_scale(fig, pos_list[0], text_list[0])
# master_length = len(text_list[0])
# print("x_unit, y units : ", x_unit, y_unit)
# annotation_df = compute_annotation_df(pos_list[1:],text_list[1:])
# print( "annotation dataframe : \n", annotation_df)
# grid = compute_grid(x_unit, y_unit)
# print("GRID :\n", grid)
# coords = find_grid_coordinates_list(pos_list, x_unit=x_unit, y_unit=y_unit)
# print("coords :", coords)
# assert len(coords) == len(pos_list)
# filled_grid = fill_grid(coords, grid)
# print("filled GRID :\n", grid)
# annotation_df["grid_coords"] = annotation_df["grid_coords"].astype('object')
# for idx in annotation_df.index :
#     grid, annotation_df = give_available_space(annotation_index= idx, annotation_df= annotation_df, grid= grid, x_unit= x_unit, y_unit=y_unit)
# annotations_obj_list = write_annotation(annotation_df,x_unit,y_unit, master_length)

# print(annotation_df)

# plt.show()