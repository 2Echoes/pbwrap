"""
This submodule contains functions to compute features related to spots detection.
"""
from bigfish.stack import check_parameter
import numpy as np
import pandas as pd
from CustomPandasFramework.PBody_project.DataFrames import newframe_Spots
from CustomPandasFramework.integrity import check_samedatashape

def compute_Spots_cell(AcquisitionId, CellId, Cell_mask, spots_dictionary : dict, cell_bbox: tuple, Pbody_label_raw) :
    """
    Parameters
    ----------
        AcquisitionId : int
            FK to acquisition table
        CellId : int
            FK to Cell table
        spots_dictionary : dict{ 'spots_type' : [(z1,y1,x1),(z2,y2,x2), ...]}
            dict object with key = spots type (such as 'rna' or 'malat1') and data is a list of 3D coordinates (z,y,x)
    """
    
    check_parameter(AcquisitionId=int, CellId=int, spots_dictionary= dict, cell_bbox= (tuple,list))
    if len(cell_bbox) != 4 : raise ValueError("Expected 4 elements in bounding box : (min_y, min_x, max_y, max_x)")

    (min_y, min_x, max_y, max_x) = cell_bbox

    types = []
    spots_coords = []
    for spot_type, coordinates_list in spots_dictionary.items() :
        if type(coordinates_list) != list : coordinates_list = list(coordinates_list)
        types.extend([spot_type] * len(coordinates_list))
        spots_coords.extend(coordinates_list)

    nbre_spots = len(spots_coords)
    if nbre_spots == 0 : return newframe_Spots()
    ids = np.arange(nbre_spots)
    Pbody_label = Pbody_label_raw.copy()
    Pbody_label = Pbody_label[min_y:max_y,min_x:max_x]
    Pbody_label[~Cell_mask] = 0

    Z,Y,X,*_ = zip(*spots_coords)
    Z,Y,X = np.array(Z), np.array(Y), np.array(X)

    spots_coords = tuple(zip(Z,Y + min_y,X + min_x))
    dataframe_ref = newframe_Spots()

    if Pbody_label.ndim == 2 : 
        Pbody_labels = Pbody_label[Y,X]
    elif Pbody_label.ndim == 3 :
        Pbody_labels = Pbody_label[Z,Y,X]
    else : raise ValueError("Pbody label has an unsupported dimension : {0}. Only 2D or 3D arrays are accepted.".format(Pbody_label.ndim))

    spots_dataframe = pd.DataFrame({
        'id' : ids,
        'AcquisitionId' : [AcquisitionId]*nbre_spots,
        'CellId' : [CellId]*nbre_spots,
        'PbodyId' : np.nan,
        'spots_coords' : spots_coords,
        'spots_type' : types,
        'Pbody_label' : Pbody_labels
    })

    check_samedatashape(dataframe_ref, spots_dataframe)
    return spots_dataframe

def compute_Spots_global(AcquisitionId, Cell_label, Nucleus_mask, Pbody_label, spots_dictionary) :
    
    types = []
    spots_coords = []
    for spot_type, coordinates_list in spots_dictionary.items() :
        if len(coordinates_list) == 0 : continue
        elif len(coordinates_list[0]) == 0 : continue
        if type(coordinates_list) != list : coordinates_list = list(coordinates_list)
        types.extend([spot_type] * len(coordinates_list))
        spots_coords.extend(coordinates_list)
    
    nbre_spots = len(spots_coords)
    if nbre_spots == 0 : return newframe_Spots()
    ids = np.arange(nbre_spots)

    Z,Y,X,*_ = zip(*spots_coords)
    Z,Y,X = np.array(Z), np.array(Y), np.array(X)

    InNucleus = Nucleus_mask[Y,X].astype(bool)

    dataframe_ref = newframe_Spots()

    if Pbody_label.ndim == 2 : 
        Pbody_labels = Pbody_label[Y,X]
    elif Pbody_label.ndim == 3 :
        Pbody_labels = Pbody_label[Z,Y,X]
    else : raise ValueError("Pbody label has an unsupported dimension : {0}. Only 2D or 3D arrays are accepted.".format(Pbody_label.ndim))

    if Cell_label.ndim == 2 : 
        Cell_labels = Cell_label[Y,X]
    elif Pbody_label.ndim == 3 :
        Cell_labels = Cell_label[Z,Y,X]
    else : raise ValueError("Pbody label has an unsupported dimension : {0}. Only 2D or 3D arrays are accepted.".format(Pbody_label.ndim))

    spots_dataframe = pd.DataFrame({
        'id' : ids,
        'AcquisitionId' : [AcquisitionId] * nbre_spots,
        'CellId' : np.nan,
        'PbodyId' : np.nan,
        'spots_coords' : spots_coords,
        'spots_type' : types,
        'cell_label' : Cell_labels,
        'Pbody_label' : Pbody_labels,
        'InNucleus' : InNucleus
    })
    filter_idx = spots_dataframe.query("cell_label != 0 or Pbody_label != 0").index

    check_samedatashape(dataframe_ref, spots_dataframe)
    return spots_dataframe.loc[filter_idx,:]