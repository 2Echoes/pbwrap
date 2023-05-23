"""
This submodule contains functions to compute features no matter the data layer.
"""

import numpy as np
from bigfish.stack import check_parameter, check_array
from .utils import unzip


def compute_signalmetrics(signal:np.ndarray, mask: np.ndarray) :
    """Compute 'min', 'max', '1 percentile', '9 percentile', 'mean', 'std' and 'median' value from signal ignoring pixels not in mask.
        
        Parameters
        ----------
            signal : np.ndarray
            mask : np.ndarray

        Returns
        -------
            signalmetrics : dict
               'min', 'max', '1 percentile', '99 percentile', 'mean' and 'median'
    """
    if mask.dtype != bool : raise TypeError("'Mask' parameter should be a ndarray with dtype = bool")

    flat:np.ndarray = signal[mask]
    signalmetrics = {
        "min" : flat.min(),
        "max" : flat.max(),
        "1 percentile" : np.percentile(flat, 1),
        "99 percentile" : np.percentile(flat, 99),
        "mean" : flat.mean(),
        "std" : flat.std(),
        "median" : np.median(flat) 
    }
    return signalmetrics


def count_spots_in_mask(spots, mask) :
    """
    Parameters
    ----------

    """

    check_parameter(spots = (list, np.ndarray), mask = np.ndarray)
    if len(spots) == 0 : return 0
    dim = len(spots[0])

    if dim == 1 : raise Exception("1D spots are not supported")

    
    if dim == 2 :
        line, col = unzip(spots)
        spot_array = np.zeros_like(mask)
        spot_array[line,col] = 1
    else : 
        plane,line,col,*_ = unzip(spots)
        #if segmentation was performed in 2D we extend constantly z-wise.
        if mask.ndim == 2 :
            mask = np.stack([mask]* (np.array(plane).max()+1))
            spot_array = np.zeros_like(mask)
        spot_array[plane,line,col] = 1

    count_array = np.logical_and(spot_array, mask)
    count = count_array.sum()

    return count



def compute_mask_area(mask: np.ndarray, unit: str = 'px', voxel_size: tuple= None)-> float:
    """
    Return the area of pbody within cell. 
    
    Parameters
    ----------
        mask: np.ndarray
            mask computed for current cell.
        unit : str
            Unit parameter can be set to either 'px' or 'nm'. If nm is selected voxel_size (z,y,x) or (y,x) has to be given as well.
        voxel_size : tuple
    """
    #Mask must be boolean
    if mask.dtype != bool : raise TypeError("'mask' parameter should be a ndarray with boolean dtype.")

    #Set pixel or nanometers
    if unit.upper() in ["PX", "PIXEL"] : return_pixel = True
    elif unit.upper() in ["NM", "NANOMETER"] and type(voxel_size) in [tuple, list] : return_pixel = False
    elif unit.upper() in ["NM", "NANOMETER"] and voxel_size == None : raise TypeError("When scale is set to nanometer, voxel_size has to be given as a tuple or a list.")
    else : raise ValueError("unit parameter incorrect should be either 'px' or 'nm'. {0} was given".format(unit))
    
    if mask.ndim != 2 : raise ValueError("Only 2D masks are supported")

    if len(voxel_size) == 2 :
        y_dim = voxel_size[0]
        x_dim = voxel_size[1]
    elif len(voxel_size) == 3 :
        y_dim = voxel_size[1]
        x_dim = voxel_size[2]
    else : raise ValueError("Inapropriate voxel_size length. Should be either 2 or 3, it is {0}".format(len(voxel_size)))


    pixel_number = mask.sum()

    if return_pixel : res = pixel_number
    else : res = pixel_number * y_dim * x_dim

    return res