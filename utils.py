import time
import numpy as np
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter, mean_filter

def from_label_get_centeroidscoords(label):
    """Returns dict{"label", "centroid"}"""

    check_parameter(label = (np.ndarray))
    centroid = regionprops_table(label, properties= ["label","centroid"])
    return centroid

def get_elmtindex(elmt,List) :
    """Returns index (position) of elmt in list
    
    Parameters
    ----------
        elmt : any
        list : list
    Returns
    -------
        res : int
    """

    check_parameter(List = (list))

    for idx in range(0,len(List)) :
        if List[idx] == elmt : yield idx


def rm_value_from_list(value,lis : list) :
    while True :
        try : lis.remove(value)
        except Exception : break
    return lis


def is_contained(list1, list2) : 
    """ Returns True if all list1 elements are in list2
    
    Parameters
    ----------
        list1 : list
        list2 : list
        
    Returns
    -------
        res : bool
        
    """

    check_parameter(list1 = (list), list2 = (list))
    truth = []

    for elmt in list1 : truth += [elmt in list2]
    res = all(truth)

    return res


def show_process_time(func):
    def wrapper(text) :
        def inner(*args, **kargs) :
            clock = time.process_time()
            res =func(*args, **kargs)
            print(text + "{0}".format(time.process_time() - clock))
            return res
        return inner
    return wrapper


def nanometer_to_pixel(value, scale) :
    if isinstance(scale, (float,int)) : scale = [scale]
    if isinstance(value, (float,int)) : value = [value]*len(scale)
    if len(value) != len(scale) : raise ValueError("value and scale must have the same dimensionality")

    return list(np.array(value) / np.array(scale))

def compute_anisotropy_coef(voxel_size) :
    """
    voxel_size : tuple (z,y,x).
    """

    if not isinstance(voxel_size, (tuple, list)) : raise TypeError("Expected voxel_size tuple or list")
    if len(voxel_size) == 2 : is_3D = False
    elif len(voxel_size) == 3 : is_3D = True
    else : raise ValueError("Expected 2D or 3D voxel, {0} element(s) found".format(len(voxel_size)))

    if is_3D :
        z_anisotropy = voxel_size[0] / voxel_size [2]
        xy_anisotropy = voxel_size[1] / voxel_size [2]
        return (z_anisotropy, xy_anisotropy, 1)

    else :
        return (voxel_size[0] / voxel_size[1], 1)
    

def inv_FWHM(FWHM) :
    """
    Returns gaussian variance for gaussian defined with Full Width at Half Maximum.
    """
    return 1/(2*np.sqrt(2* np.log(2))) * FWHM

def inv_FWTM(FWTM) :
    """
    Returns gaussian variance for gaussian defined with Full Width at Tenth Maximum.
    """
    return 1/(2*np.sqrt(2* np.log(2))) * FWTM


def gaussian_kernel_size(object_size_nm, voxel_size, width= 'FWHM') :
    """
    Computes the kernel size of a Gaussian function so that either Full Width at Half Maximum (width = 'FWHM') or Full Width at Tenth Maximum (width = 'FWTM') corresponds to the object dimension.
    Object_size_nm should be coherent in type (and length if tuples or lists are passed).
    Always returns a list with object dim number of element.
    """
    
    if width == 'FWHM' : variance_func = inv_FWHM
    elif width == 'FWTM' : variance_func = inv_FWTM
    else : raise ValueError("with should either be 'FWHM' or 'FWTM'. It is {0}".format(width))

    if isinstance(object_size_nm, (tuple, list)) and isinstance(voxel_size, (tuple, list)):
        if len(object_size_nm) != len(voxel_size) : raise ValueError("Length of object_size_nm and voxel_size parameters should be coherent.")
        return [variance_func(obj_size_nm / scale_factor) for obj_size_nm, scale_factor in zip(object_size_nm, voxel_size)]
    elif isinstance(object_size_nm, (float,int)) and isinstance(voxel_size, (tuple, list)): 
        return [variance_func(object_size_nm / voxel_size)]
    else : raise TypeError("object size and voxel_size parameters should be tuple or float like object and be coherent.")

def compute_voxel_volume(voxel) : 
    return np.array(voxel).prod()


def _get_varname(variable):
    for name in globals():
        if id(globals()[name]) == id(variable):
            return name
    for name in locals():
        if id(locals()[name]) == id(variable):
            return name
    return None


def get_centroids_list(clusters_df) :

    """
    clusters_list should be a pd.DataFrame with ['z', 'y', 'x'] or ['y', 'x'] keys.
    """

    if 'y' in clusters_df.columns and 'x' in clusters_df.columns :
        if 'z' in clusters_df.columns : keys = [clusters_df['z'], clusters_df['y'], clusters_df['x']]
        else : keys = [clusters_df['y'], clusters_df['x']]
    else : raise ValueError("Expected keys : ['z', 'y', 'x'] or ['y', 'x']")

    return list(zip(*keys))

def get_centroids_array(cluster_df) :

    return np.array(get_centroids_list(cluster_df), dtype= int)