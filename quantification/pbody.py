import numpy as np
import pandas as pd
import time
import CustomPandasFramework.PBody_project.DataFrames as Dataframe
from CustomPandasFramework.integrity import check_samedatashape
from .measures import count_spots_in_mask, count_spots_in_masks_list, count_rna_close_pbody, count_rna_close_pbody_list, count_rna_close_pbody_global
from skimage.segmentation import find_boundaries
from skimage.measure import regionprops_table
from pbwrap.integrity import check_parameter

def compute_Pbody(AcquisitionId: int, Pbody_label: np.ndarray, cell_label: np.ndarray, rna_coords, malat1_coords) :
    """
    Compute Pbody DF during analysis pipeline.
    Note : It is important that the Pbody_label given is the same as the one used during cells computation (fov).
    """
    Pbody_dim = Pbody_label.ndim
    cell_dim = cell_label.ndim
    Pbody_dictionary = compute_Pbody_dictionary(Pbody_label, rna_coords, malat1_coords)
    nbre_pbody = len(Pbody_dictionary["label"])
    ids = np.arange(nbre_pbody)
    AcquisitionIds = [AcquisitionId] * nbre_pbody
    CellIds = np.nan # Computed during 'CustomPandasFramework.Pbody_project.Pbody_AddCellFK' call
    centroids_coordinates = list(zip(*(np.array(Pbody_dictionary["centroid-{0}".format(n)]).round().astype(int) for n in range(0,Pbody_dim)))) 

    if Pbody_dim == 2 :
        areas  = Pbody_dictionary["area"]
        volumes = np.nan
    elif Pbody_dim == 3 :
        areas = len(Pbody_dictionary["boundary"])
        volumes = Pbody_dictionary["area"] #TODO : test that this has been implemented as stated on the internet.    

    Y,X = np.array(Pbody_dictionary["centroid-0"]).round().astype(int), np.array(Pbody_dictionary["centroid-1"]).round().astype(int)
    
    if cell_dim == 2 : cell_labels = cell_label[Y,X]
    else : raise ValueError("Only 2D arrays are supported for Cell label")
    DF_1000nm = pd.DataFrame(columns= ["label", "rna_closer_1000nm"], data=list(zip(*Pbody_dictionary["pbody_closer_than_1000_nm"]))) 
    DF_1500nm = pd.DataFrame(columns= ["label", "rna_closer_1500nm"], data=list(zip(*Pbody_dictionary["pbody_closer_than_1500_nm"]))) 
    DF_2000nm = pd.DataFrame(columns= ["label", "rna_closer_2000nm"], data=list(zip(*Pbody_dictionary["pbody_closer_than_2000_nm"])))

    DF = pd.merge(DF_2000nm, DF_1500nm,how= 'left', on= 'label') 
    DF = pd.merge(DF, DF_1000nm,how= 'left', on= 'label') 

    res_DataFrame = pd.DataFrame({
        "id" : ids,
        "AcquisitionId" : AcquisitionIds,
        "CellId" : CellIds,
        "centroid_coordinates" : centroids_coordinates,
        "area" : areas,
        "volume" : volumes,
        "label" : Pbody_dictionary["label"],
        "cell_label" : cell_labels,
        "rna_count" : Pbody_dictionary["rna_count"],
        "malat1_count" : Pbody_dictionary["malat1_count"]
    })
    res_DataFrame = pd.merge(res_DataFrame, DF, 'left', on= 'label')
    datashape_ref = Dataframe.newframe_Pbody()
    check_samedatashape(res_DataFrame, datashape_ref)
    res_DataFrame = res_DataFrame.query('cell_label != 0')

    print(Pbody_label.sum())
    return res_DataFrame



def compute_Pbody_dictionary(Pbody_label: np.ndarray, rna_coords, malat_coords, voxelsize = (300,103,103)):
    """
    From Pbody_label (fov) computes features and return dict object. 
    Each item is a list of feature where each element is the feature value for 1 region in the label.

    Keys
    ----
        'label' : int
            label of the region
        'centroid-x' : float
            coordinate of the region centroid on the x axis. (Ex : for 2D im centroid-0, centroid-1 are computed...)
        'area' : int
            area of the mask in pixel
        'boundary_coordinates' : list[tuple]
            Each element contains a tuple defining the border of one pbody. Each of these tuples contains the coordinates of every point in the boundary of the pbody.
    """

    if not isinstance(Pbody_label, np.ndarray) : raise TypeError('Pbody_label should be of ndarray type. It is {0}'.format(type(Pbody_label)))
    Pbody_dictionary = regionprops_table(Pbody_label, properties= ['label', 'centroid','area','bbox'])
    masks_list = [Pbody_label == label for label in Pbody_dictionary['label']]
    rna_count = count_spots_in_masks_list(rna_coords, masks_list)
    malat1_count = count_spots_in_masks_list(malat_coords, masks_list)

    # clock = time.process_time()
    # pbody_closer_than_1000_nm = [count_rna_close_pbody(pbody_mask= pbody_mask, spots_coords= rna_coords, distance_nm= 1000, voxel_size= voxelsize) for pbody_mask in masks_list]
    # print("time for 1st individual computation : {0}".format(time.process_time() - clock))
    clock = time.process_time()
    # pbody_closer_than_1500_nm = [count_rna_close_pbody(pbody_mask= pbody_mask, spots_coords= rna_coords, distance_nm= 1500, voxel_size= voxelsize) for pbody_mask in masks_list]
    # pbody_closer_than_2000_nm = [count_rna_close_pbody(pbody_mask= pbody_mask, spots_coords= rna_coords, distance_nm= 2000, voxel_size= voxelsize) for pbody_mask in masks_list]
    Pbody_dictionary["rna_count"] = rna_count
    Pbody_dictionary["malat1_count"] = malat1_count
    Pbody_dictionary["pbody_closer_than_1000_nm"] = count_rna_close_pbody_global(pbody_label= Pbody_label, spots_coords= rna_coords, distance_nm= 1000, voxel_size= voxelsize)
    Pbody_dictionary["pbody_closer_than_1500_nm"] = count_rna_close_pbody_global(pbody_label= Pbody_label, spots_coords= rna_coords, distance_nm= 1500, voxel_size= voxelsize)
    Pbody_dictionary["pbody_closer_than_2000_nm"] = count_rna_close_pbody_global(pbody_label= Pbody_label, spots_coords= rna_coords, distance_nm= 2000, voxel_size= voxelsize)
    print("time for 3 global computation : {0}".format(time.process_time() - clock))

    return Pbody_dictionary


def boundary_coordinates(mask: np.ndarray, to_tuple = True) -> tuple :

    """
    Is not used due to a too long computational time : ~ 2mins per fov
    from mask returns a tuple of coordinates where each element are the coordinates of the mask boundaries.
    If to_tuple = False : returns list instead.
    """

    boundary = find_boundaries(mask, background= 0, mode = 'inner')
    boundary_coordinates = np.argwhere(boundary).tolist()

    if to_tuple: return tuple(map(tuple, boundary_coordinates))
    else : return boundary_coordinates

