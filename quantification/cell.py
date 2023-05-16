import numpy as np
import warnings
from scipy.ndimage import distance_transform_edt
from bigfish.stack import mean_projection, maximum_projection, check_parameter
from ..utils import rm_value_from_list
from .utils import unzip


def nucleus_signal_metrics(cell, channel, projtype = 'mip') :
    """Returns dict containing signal related measures : 'min', 'max', '1 percentile', '9 percentile', 'mean' and 'median'.
      Computed from channel signal in cell's nucleus mask. Signal measures are computed from 2D cell, so channel is projected z-wise according to projtype (provided channel is 3D).
    
        Parameters
        ----------
            cell : dict
                Dictionary computed from bigFish.multistack.extract_cell
            channel : np.ndarray
                Channel from which intensity is computed
            projtype : str
                can either be 'mip' or 'mean'.

        Returns
        -------
            mean_sig : float
        
    """
    min_y, min_x, max_y, max_x = cell["bbox"]
    channel_cropped = channel[:, min_y:max_y, min_x:max_x]
 

    if channel.ndim == 3 :
        if projtype == 'mip' : 
            channel_crop_proj = maximum_projection(channel_cropped)
        elif projtype == 'mean' :
            channel_crop_proj = mean_projection(channel_cropped)
    
    nucleus_mask = cell["nuc_mask"]
    metrics = compute_signalmetrics(channel_crop_proj, nucleus_mask)
    return metrics


def OutOfNucleus_signal_metrics(cell, channel, projtype = 'mip') :
    """Returns mean signal from channel in computed cell. Mean signal is computed from 2D cell, so channel is projected z-wise according to projtype provided channel is 3D.
    
        Parameters
        ----------
            cell : dict
                Dictionary computed from bigFish.multistack.extract_cell
            channel : np.ndarray
                Channel from which intensity is computed
            projtype : str
                can either be 'mip' or 'mean'.

        Returns
        -------
            mean_sig : float
        
    """
    min_y, min_x, max_y, max_x = cell["bbox"]
    channel_cropped = channel[:, min_y:max_y, min_x:max_x]
 

    if channel.ndim == 3 :
        if projtype == 'mip' : 
            channel_crop_proj = maximum_projection(channel_cropped)
        elif projtype == 'mean' :
            channel_crop_proj = mean_projection(channel_cropped)
    
    OutOfNucleus_mask = np.logical_not(cell["nuc_mask"])
    metrics = compute_signalmetrics(channel_crop_proj, OutOfNucleus_mask)
    return metrics



def compute_signalmetrics(signal, mask) :
    """Compute 'min', 'max', '1 percentile', '9 percentile', 'mean' and 'median' value from signal ignoring pixels not in mask.
        
        Parameters
        ----------
            signal : np.ndarray
            mask : np.ndarray

        Returns
        -------
            signalmetrics : dict
               'min', 'max', '1 percentile', '99 percentile', 'mean' and 'median'
    """

    signal[mask != 1] = -999
    flat = list(signal.flatten())
    flat = np.array(rm_value_from_list(-999, flat))
    signalmetrics = {
        "min" : flat.min(),
        "max" : flat.max(),
        "1 percentile" : np.percentile(flat, 1),
        "99 percentile" : np.percentile(flat, 99),
        "mean" : flat.mean(),
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



def compute_mask_area(pbody_mask: np.ndarray, unit: str = 'px', voxel_size: tuple= None)-> float:
    """
    Return the area of pbody within cell. 
    
    Parameters
    ----------
        pbody_mask: np.ndarray
            mask computed for current cell.
        unit : str
            Unit parameter can be set to either 'px' or 'nm'. If nm is selected voxel_size (z,y,x) or (y,x) has to be given as well.
        voxel_size : tuple
    """
    if unit.upper() in ["PX", "PIXEL"] : return_pixel = True
    elif unit.upper() in ["NM", "NANOMETER"] and type(voxel_size) in [tuple, list] : return_pixel = False
    elif unit.upper() in ["NM", "NANOMETER"] and voxel_size == None : raise TypeError("When scale is set to nanometer, voxel_size has to be given as a tuple or a list.")
    else : raise ValueError("unit parameter incorrect should be either 'px' or 'nm'. {0} was given".format(unit))
    if pbody_mask.ndim != 2 : raise ValueError("Only 2D masks are supported")

    if len(voxel_size) == 2 :
        y_dim = voxel_size[0]
        x_dim = voxel_size[1]
    elif len(voxel_size) == 3 :
        y_dim = voxel_size[1]
        x_dim = voxel_size[2]
    else : raise ValueError("Inapropriate voxel_size length. Should be either 2 or 3, it is {0}".format(len(voxel_size)))

    _, count = np.unique(pbody_mask, return_counts= True)
    del _
    assert len(count) == 2
    pixel_number = count[1]

    if return_pixel : res = pixel_number
    else : res = pixel_number * y_dim * x_dim

    return res



def count_rna_close_pbody(pbody_mask: np.ndarray, spots_coords: 'list[tuple]', distance_nm: float, voxel_size: 'tuple[float]')-> int :
    """
    Count number of RNA (spots) closer than 'distance_nm' from a p-body (mask).
    """
    
    check_parameter(pbody_mask = (np.ndarray), spots_coords = (list, np.ndarray), distance_nm = (int, float), voxel_size = (tuple, list))

    if pbody_mask.ndim != 2: raise ValueError("Unsupported p_body mask dimension. Only 2D arrays are supported.")
    if type(spots_coords) == np.ndarray : spots_coords = list(spots_coords)
    if len(voxel_size) == 3 :
        y_scale = voxel_size[1]
        x_scale = voxel_size[2]
    elif len(voxel_size) == 2 :
        y_scale = voxel_size[0]
        x_scale = voxel_size[1]
    else : raise ValueError("Incorrect voxel_size length should be either 2 or 3. {0} was given".format(len(voxel_size)))


    frompbody_distance_map = distance_transform_edt(np.logical_not(pbody_mask), sampling= [y_scale, x_scale])
    rna_distance_map = np.ones_like(pbody_mask) * -999
    if len(spots_coords[0]) == 2 :
        y_coords, x_coords = unzip(spots_coords)
    elif len(spots_coords[0]) == 3 :
        z_coords, y_coords, x_coords = unzip(spots_coords)
        del z_coords
    else : 
        z_coords, y_coords, x_coords,*_ = unzip(spots_coords)
        del z_coords,_
    rna_distance_map[y_coords, x_coords] = frompbody_distance_map[y_coords, x_coords] # This distance maps gives the distance of each RNA to the closest p-body
    count_map = rna_distance_map[rna_distance_map >= 0] <= distance_nm
    
    values,count = np.unique(count_map, return_counts= True)
    if not True in values : 
        count = 0
    else:
        index = list(values).index(True)
        count = count[index]
    
    return count

        