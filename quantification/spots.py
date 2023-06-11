"""
This submodule contains functions to compute features related to spots detection.
"""
from bigfish.stack import check_parameter
import numpy as np
import pandas as pd
from CustomPandasFramework.PBody_project.DataFrames import newframe_Spots
from CustomPandasFramework.operations import check_samedatashape

def compute_Spots(AcquisitionId, CellId, spots_dictionary : dict, cell_bbox: tuple) :
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
    if len(cell_bbox) != 4 : raise ValueError("Expected 4 elements en bounding box : (min_y, min_x, max_y, max_x)")

    (min_y, min_x, max_y, max_x) = cell_bbox

    nbre_spots = len(spots_coords.values())
    ids = np.arange(nbre_spots)

    types = []
    spots_coords = []
    for spot_type, coordinates_list in spots_dictionary.items() :
        types += [spot_type] * len(coordinates_list)
        spots_coords += coordinates_list


    spots_coords = spots_dictionary.values()
    Z,Y,X = zip(*spots_coords)
    Z,Y,X = np.array(Z), np.array(Y), np.array(X)
    Y += min_y
    X +=  min_x 
    spots_coords = tuple(zip(Z,Y,X))
    types = spots_dictionary.keys()
    dataframe_ref = newframe_Spots()

    spots_dataframe = pd.DataFrame({
        'AcquisitionId' : [AcquisitionId]*nbre_spots,
        'CellId' : [CellId]*nbre_spots,
        'id' : ids,
        'spots_coords' : spots_coords,
        'type' : types
    })

    check_samedatashape(dataframe_ref, spots_dataframe)

    return spots_dataframe
