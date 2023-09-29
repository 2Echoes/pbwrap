import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.DataFrames as Dataframe
from CustomPandasFramework.integrity import check_samedatashape
from .measures import count_spots_in_mask, count_spots_in_masks_list, count_rna_close_pbody, count_rna_close_pbody_list, count_rna_close_pbody_global
from skimage.segmentation import find_boundaries
from skimage.measure import regionprops_table
from pbwrap.integrity import check_parameter

def compute_Pbody(AcquisitionId: int, Pbody_label: np.ndarray, cell_label: np.ndarray, nucleus_mask: np.ndarray, rna_coords, malat1_coords, distance = [0, 100,200,400,600,800,1000,1500,2000]) :
    """
    Compute Pbody DF during analysis pipeline.
    Note : It is important that the Pbody_label given is the same as the one used during cells computation (fov).
    """
    Pbody_dim = Pbody_label.ndim
    cell_dim = cell_label.ndim
    Pbody_dictionary = compute_Pbody_dictionary(Pbody_label, rna_coords, malat1_coords, distance=distance)
    nbre_pbody = len(Pbody_dictionary["label"])
    ids = np.arange(nbre_pbody)
    AcquisitionIds = [AcquisitionId] * nbre_pbody
    CellIds = np.nan # Computed during 'CustomPandasFramework.Pbody_project.Pbody_AddCellFK' call
    centroids_coordinates = list(zip(*(np.array(Pbody_dictionary["centroid-{0}".format(n)]).round().astype(int) for n in range(0,Pbody_dim)))) 
    InNucleus = nucleus_mask[Pbody_dictionary["centroid-0"].round().astype(int), Pbody_dictionary["centroid-1"].round().astype(int)]
    if Pbody_dim == 2 :
        areas  = Pbody_dictionary["area"]
        volumes = np.nan
    elif Pbody_dim == 3 :
        areas = len(Pbody_dictionary["boundary"])
        volumes = Pbody_dictionary["area"] #TODO : test that this has been implemented as stated on the internet.    

    Y,X = np.array(Pbody_dictionary["centroid-0"]).round().astype(int), np.array(Pbody_dictionary["centroid-1"]).round().astype(int)
    
    if cell_dim == 2 : cell_labels = cell_label[Y,X]
    else : raise ValueError("Only 2D arrays are supported for Cell label")
    DF = pd.DataFrame({'label' : Pbody_dictionary['label']})
    for dist in distance :
        DF = pd.merge(DF, pd.DataFrame(columns= ['label', 'rna {0}nm count'.format(dist)], data= zip(Pbody_dictionary['rna {0} nm'.format(dist)].index, Pbody_dictionary['rna {0} nm'.format(dist)])), how= 'left', on= 'label')
        DF = pd.merge(DF, pd.DataFrame(columns= ['label', 'malat1 {0}nm count'.format(dist)], data= zip(Pbody_dictionary['malat1 {0} nm'.format(dist)].index, Pbody_dictionary['malat1 {0} nm'.format(dist)])), how= 'left', on= 'label')
    DF = DF.fillna(0)

    res_DataFrame = pd.DataFrame({
        "id" : ids,
        "AcquisitionId" : AcquisitionIds,
        "CellId" : CellIds,
        "centroid_coordinates" : centroids_coordinates,
        "area" : areas,
        "volume" : volumes,
        "label" : Pbody_dictionary["label"],
        "cell_label" : cell_labels,
        "InNucleus" : InNucleus
    })

    res_DataFrame = pd.merge(res_DataFrame, DF, 'left', on= 'label')
    res_DataFrame = res_DataFrame.query('cell_label != 0')
    datashape_ref = Dataframe.newframe_Pbody(distance=distance)
    check_samedatashape(res_DataFrame, datashape_ref)
    return res_DataFrame



def compute_Pbody_dictionary(Pbody_label: np.ndarray, rna_coords, malat_coords, distance, voxelsize = (300,103,103)):
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

    rna_count_dictionary = count_rna_close_pbody_global(pbody_label= Pbody_label, spots_coords= rna_coords, distance_nm= distance, voxel_size= voxelsize, spot_type= 'rna')
    Pbody_dictionary.update(rna_count_dictionary)
    malat1_count_dictionary = count_rna_close_pbody_global(pbody_label= Pbody_label, spots_coords= malat_coords, distance_nm= distance, voxel_size= voxelsize, spot_type= 'malat1')
    Pbody_dictionary.update(malat1_count_dictionary)


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

