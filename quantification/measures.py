"""
This submodule contains functions to compute features no matter the data layer.
"""

import numpy as np
import pandas as pd
import bigfish.stack as stack
from bigfish.stack import check_parameter, check_array
from .utils import unzip
from scipy.ndimage import distance_transform_edt
from scipy.signal import fftconvolve
from ..utils import check_parameter
from ..utils import nanometer_to_pixel as utils_nanometer_to_pixel

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
        "median" : np.median(flat),
        "sum" : np.sum(flat) 
    }
    return signalmetrics


def count_spots_in_mask(spots, mask: np.ndarray)->int :
    """
    Parameters
    ----------

    """

    check_parameter(spots = (list, np.ndarray), mask = np.ndarray)
    if len(spots) == 0 : return 0
    dim = len(spots[0])
    if mask.dtype != bool : mask = mask.astype(bool)


    
    if dim == 2 :
        line, col = unzip(spots)
        count = mask[line, col].sum()
    elif dim == 1 : 
        raise Exception("1D spots are not supported")

    else : 
        plane,line,col,*_ = unzip(spots)
        #if segmentation was performed in 2D we extend constantly z-wise.
        if mask.ndim == 2 :
            count = mask[line,col].sum()
        else : count = mask[plane,line, col].sum()


    return count

def count_spots_in_masks_list(spots, masks: 'list[np.ndarray]')-> 'list[int]' :
    """"
    Similar to count_spots_in_mask but with the a list of masks. Same spots coords are used to count in every masks. But way faster than looping over a masks list while indexing count_spots_in_masks results
    Returns a list of counts of length len(masks)
    """

    if type(masks) == list :
        masks = np.array(masks)
    elif type(masks) == np.ndarray :
        if not masks.ndim in [3,4] : raise ValueError("Unsupported dimension for masks parameter, ndim should be 3 for 2D maskS and 4 for 3D maskS. It is {0}".format(masks.ndim))
    else : raise TypeError("Masks should be a list or np.ndarray. It is {0}".format(type(masks)))
    
    if len(spots) == 0 : return [0]*masks.shape[0]
    dim = len(spots[0])

    if dim == 2 :
        line, col = unzip(spots)
        count = masks[:,line, col].sum(axis= 1)
    elif dim == 1 : 
        raise Exception("1D spots are not supported")

    else : 
        plane,line,col,*_ = unzip(spots)
        #if segmentation was performed in 2D we extend constantly z-wise.
        if masks.ndim == 3 :
            count = masks[:, line,col].sum(axis= 1)
        else : count = masks[:,plane,line, col].sum(axis= 1)

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

    if not return_pixel :
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

def nucleus_signal_metrics(cell, channel, projtype = 'mip', use_cell_mask= False) :
    """
    Returns dict containing signal related measures : 'min', 'max', '1 percentile', '9 percentile', 'mean' and 'median'.
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
            channel_crop_proj = stack.maximum_projection(channel_cropped)
        elif projtype == 'mean' :
            channel_crop_proj = stack.mean_projection(channel_cropped)
    
    if use_cell_mask : nucleus_mask = cell["cell_mask"]
    else : nucleus_mask = cell["nuc_mask"]

    metrics = compute_signalmetrics(channel_crop_proj, nucleus_mask)
    return metrics



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
    if len(spots_coords) == 0 : return 0
    if len(spots_coords[0]) == 2 :
        y_coords, x_coords = unzip(spots_coords)
    elif len(spots_coords[0]) == 3 :
        z_coords, y_coords, x_coords = unzip(spots_coords)
        del z_coords
    else : 
        z_coords, y_coords, x_coords,*_ = unzip(spots_coords)
        del z_coords,_
    rna_distance_map[y_coords, x_coords] = frompbody_distance_map[y_coords, x_coords] # This distance maps gives the distance of each RNA to the closest p-body
    # count_map = rna_distance_map[rna_distance_map >= 0] <= distance_nm
    count_map = np.logical_and(rna_distance_map >= 0, rna_distance_map <= distance_nm)
    count = np.sum(count_map, axis= 1)

    return count


def count_rna_close_pbody_list(list_pbody_mask: np.ndarray, spots_coords: 'list[tuple]', distance_nm: float, voxel_size: 'tuple[float]')-> 'list[int]' :
    
    if len(voxel_size) == 3 :
        y_scale = voxel_size[1]
        x_scale = voxel_size[2]
    elif len(voxel_size) == 2 :
        y_scale = voxel_size[0]
        x_scale = voxel_size[1]
    else : raise ValueError("Incorrect voxel_size length should be either 2 or 3. {0} was given".format(len(voxel_size)))
    pbody_masks = np.array(list_pbody_mask)
    # frompbody_distance_map = distance_transform_edt(np.logical_not(pbody_masks), sampling= [0, y_scale, x_scale]) # 1e15 used to make the computation unrelated to the z axis which has no meaning here.
    frompbody_distance_map = np.array([distance_transform_edt(np.logical_not(pbody_mask), sampling= [y_scale, x_scale]) for pbody_mask in pbody_masks])
    rna_distance_map = np.ones_like(pbody_masks) * -999
    
    if len(spots_coords) == 0 : return 0
    if len(spots_coords[0]) == 2 :
        y_coords, x_coords = unzip(spots_coords)
    elif len(spots_coords[0]) == 3 :
        z_coords, y_coords, x_coords = unzip(spots_coords)
        del z_coords
    else : 
        z_coords, y_coords, x_coords,*_ = unzip(spots_coords)
        del z_coords,_
    z_coords = np.arange(len(list_pbody_mask))
    rna_distance_map[:, y_coords, x_coords] = frompbody_distance_map[:, y_coords, x_coords] # This distance maps gives the distance of each RNA to the closest p-body
    count_map = np.logical_and(rna_distance_map >= 0, rna_distance_map <= distance_nm)
    counts = np.sum(count_map, axis= (2,1))
    return counts

def count_rna_close_pbody_global(pbody_label: np.ndarray, spots_coords: 'list[tuple]', distance_nm: float, voxel_size: 'tuple[float]', spot_type:str ='spot')-> dict :
    """
    Count number of RNA (spots) closer than 'distance_nm' from a p-body (mask). 
    Distance_nm argument can also be a list[float like] in such a case one measure will be computed for each element of the list. Therefor output will be a list of length = len(distance_nm)
    
    Returns
    -------
        count_dictionary : dictionary
            dictionary with length = len(distance_nm) or 1 if distance_nm is float. Keys are element of distance_nm.
            Such as {'distance1' : count_d1, 'distance2' : count_d2....}
    """
    
    check_parameter(pbody_label = (np.ndarray), spots_coords = (list, np.ndarray), distance_nm = (int, float, list), voxel_size = (tuple, list))

    if pbody_label.ndim != 2: raise ValueError("Unsupported p_body mask dimension. Only 2D arrays are supported.")
    if type(spots_coords) == np.ndarray : spots_coords = list(spots_coords)
    if len(voxel_size) == 3 :
        y_scale = voxel_size[1]
        x_scale = voxel_size[2]
    elif len(voxel_size) == 2 :
        y_scale = voxel_size[0]
        x_scale = voxel_size[1]
    else : raise ValueError("Incorrect voxel_size length should be either 2 or 3. {0} was given".format(len(voxel_size)))

    res = {'{0} {1} nm'.format(spot_type, distance): pd.Series(dtype= float) for distance in distance_nm}

    if len(spots_coords) == 0 :
        return res
    elif len(spots_coords[0]) == 0 :
        return res
    elif len(spots_coords[0]) == 2 :
        y_coords, x_coords = unzip(spots_coords)
    elif len(spots_coords[0]) == 3 :
        z_coords, y_coords, x_coords = unzip(spots_coords)
        del z_coords
    else : 
        z_coords, y_coords, x_coords,*_ = unzip(spots_coords)
        del z_coords,_

    #Constructing a frame to enable counting spots with same coordinates.
    spots_number_frame = pd.DataFrame({"spots_coords" : list(zip(y_coords,x_coords)), "id" : np.arange(len(y_coords))})
    spots_number_frame = spots_number_frame.groupby(["spots_coords"])["id"].count().rename("count")
    #Distance map
    shape = pbody_label.shape
    pbody_mask = pbody_label.copy().astype(bool)
    frompbody_distance_map, indices = distance_transform_edt(np.logical_not(pbody_mask), sampling= [y_scale, x_scale], return_indices= True)  #assert z_scale > y,x scale because pbody segmentation is done in 2D so we want to make sure the shortest distance to pbodies will always be found in the plane.
    rna_distance_map = np.ones(shape) * -999
    rna_distance_map[y_coords, x_coords] = frompbody_distance_map[y_coords, x_coords] # This distance maps gives the distance of each RNA to the closest p-body
    
    #Counting
    if isinstance(distance_nm, (int, float)) : distance_nm = [distance_nm]
    count_maps = [np.logical_and(rna_distance_map >= 0, rna_distance_map <= distance) for distance in distance_nm]
    res = {
        '{0} {1} nm'.format(spot_type, distance):
        pd.DataFrame(
        columns= ['spots_coords', 'label'],
        data = [((Y,X), pbody_label[indices[0,Y,X], indices[1,Y,X]]) for Y,X in zip(np.nonzero(count_map)[0],np.nonzero(count_map)[1])]
        ).set_index(['spots_coords']).join(spots_number_frame, on= 'spots_coords').groupby('label')['count'].sum()
    for count_map,distance in zip(count_maps, distance_nm)
    }

    return res



def reconstruct_boolean_signal(image_shape, spot_list: list):
    signal = np.zeros(image_shape, dtype= bool)
    if len(spot_list) == 0 : return signal
    Z, Y, X = list(zip(*spot_list))
    signal[Z,Y,X] = True

    return signal



def spots_colocalisation(image_shape, spot_list1:list, spot_list2:list, distance: float, voxel_size)-> int :
    """
    Return number of spots from spot_list1 located closer(large) than distance to at least one spot of spot_list2.

    Parameters
    ----------
        image_shape : tuple
        spot_list1 : list
        spot_list2 : list
        distance : nanometer
            distance in nanometer.
        voxel_size : (z,y,x) tuple
    """

    if len(image_shape) != 3 : raise ValueError("Image shape length should be 3 : (Z,Y,X)")
    if len(spot_list1) == 0 or len(spot_list2) == 0 : return np.NaN

    signal2 = reconstruct_boolean_signal(image_shape, spot_list2)
    mask = np.logical_not(signal2)
    distance_map = distance_transform_edt(mask, sampling= voxel_size)
    Z,Y,X = zip(*spot_list1)

    count = (distance_map[Z,Y,X] <= distance).sum()
    return count

def cluster_localisation(clusters_dataframe: pd.DataFrame) :
    return clusters_dataframe["spot_number"].sum()



###Attempt at counting number of spots within given radius of a pixel, for all pixels
def _reconstruct_spot_signal(image_shape, spot_list: list):
    """
    Create a map where each pixel value correspond to the number of spots located in this position.
    """
    signal = np.zeros(image_shape, dtype= int)
    unique_list, counts = np.unique(spot_list, return_counts= True, axis=0)
    Z, Y, X = list(zip(*unique_list))
    signal[Z,Y,X] = counts

    return signal


def _create_counting_kernel(radius_nm, voxel_size) :

    max_pixel_distance = int(max(utils_nanometer_to_pixel(radius_nm, voxel_size)))
    kernel = np.ones(shape=(2*max_pixel_distance+1 ,2*max_pixel_distance+1, 2*max_pixel_distance+1)) #always odd number so middle is always at [pixel_radius-1, pixel_radius-1, pixel_radius-1]
    kernel[max_pixel_distance, max_pixel_distance, max_pixel_distance] = 0
    kernel = distance_transform_edt(kernel, sampling= voxel_size) <= radius_nm
    
    return kernel.astype(int)


def _spot_count_map(spots_array, radius_px, voxel_size) :
    """
    Create a map where each pixel value correspond to the number of spots closer than radius to the position.
    """

    kernel = _create_counting_kernel(radius_px, voxel_size)
    map = fftconvolve(spots_array, kernel, mode= 'same')

    return np.round(map).astype(int)

    
def nanometer_to_pixel(value, scale) :
    print('depreciated : function moved to pbwrap.utils')
    if isinstance(scale, (float,int)) : scale = [scale]
    if isinstance(value, (float,int)) : value = [value]*len(scale)
    if len(value) != len(scale) : raise ValueError("value and scale must have the same dimensionality")

    return list(np.array(value) / np.array(scale))
    

def spots_multicolocalisation(spots_list, anchor_list, radius_nm, image_shape, voxel_size) :

    """
    Compute the number of spots from spots_list closer than radius to a spot from anchor_list. Each spots_list spots will be counted as many times as there are anchors close enough.
    Note that the radius in nm is converted to pixel using voxel size, and rounded to the closest int value.
    
    Example in 2D
    --------

    >>> Anchors         Spots           Radius (2px)    Count
    >>> 0 0 0 0 0 0     0 X 0 0 X 0       1             0 1 0 0 0 0
    >>> 0 X 0 0 0 0     X 0 0 X 0 0     1 1 1           1 0 0 0 0 0
    >>> 0 X 0 0 0 0     X X 0 0 0 0       1             1 2 0 0 0 0     --> 5
    >>> 0 0 0 0 X 0     0 0 X 0 0 0                     0 0 0 0 0 0
    >>> 0 0 0 0 0 0     0 0 0 X 0 0                     0 0 0 0 0 0

    Parameters
    ----------
    spots_list : list
    anchor_list : list
    radius_nm : int, float
    image_shape : tuple (Z, Y, X)
    voxel_size : tuple (Z, Y, X)
    
    Returns
    -------
    Returns the list of neighbouring spot number to 'spots_list'.
    """

    check_parameter(spots_list= (list, np.ndarray), anchor_list= (list, np.ndarray), radius_nm = (int, float), image_shape= (tuple,list), voxel_size= (tuple, list))
    if len(image_shape) != 3 : raise ValueError("Only 3D colocalisation is supported, 'image_shape' should be (Z,Y,X).")
    if len(voxel_size) != 3 : raise ValueError("Only 3D colocalisation is supported, 'voxel_size' should be (Z,Y,X).")
    if voxel_size[1] != voxel_size[2] : raise ValueError("Unsupported anisotropy in xy plan. (yscale != xscale)")
    if len(spots_list) == 0 or len(anchor_list) == 0 : return 0

    anchor_array = _reconstruct_spot_signal(image_shape=image_shape, spot_list=anchor_list)
    count_map = _spot_count_map(anchor_array, radius_px=radius_nm, voxel_size=voxel_size)
    Z,Y,X = list(zip(*spots_list))

    return list(count_map[Z,Y,X])


def closest_spot_distance(coordinates_list, spots_list, image_shape, voxel_size) :
    """
    Compute distance to closest spot (from spots_list) for each point in coordinates_list.
    """

    check_parameter(spots_list= (list, np.ndarray), coordinates_list= (list, np.ndarray), image_shape= (tuple,list), voxel_size= (tuple, list))
    if len(image_shape) != 3 : raise ValueError("Only 3D colocalisation is supported, 'image_shape' should be (Z,Y,X).")
    if len(voxel_size) != 3 : raise ValueError("Only 3D colocalisation is supported, 'voxel_size' should be (Z,Y,X).")
    if voxel_size[1] != voxel_size[2] : raise ValueError("Unsupported anisotropy in xy plan. (yscale != xscale)")
    if len(spots_list) == 0 or len(coordinates_list) == 0 : return 0

    spot_array = _reconstruct_spot_signal(image_shape=image_shape, spot_list= spots_list).astype(bool)
    distance_map = distance_transform_edt(np.logical_not(spot_array), sampling= voxel_size)
    Z,Y,X = list(zip(*coordinates_list))

    return list(distance_map[Z,Y,X])
