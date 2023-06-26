import numpy as np
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter


def unzip(lis:list)-> 'list[list]' :
    """from a list of coordinates return a list of list where each rows has all coordinate from an axis
    
    Parameters
    ----------
        lis : list, np.ndarray
        
    Returns
    -------
        res : list
    """
    if type(lis) == np.ndarray : lis = list(lis)


    res = [list(dim) for dim in zip(*lis)]
    # for dim in zip(*lis):
    #     res += [list(dim)]
    return res

def spots_z_proj(spots) :
    """Project z wize 3D coords (zyx) to 2D coords(yx). Basically it consists in removing z for coords.
    
    Parameters
    ----------
        spots : list, np.ndarray
            list of spots with coords (zyx)
    """

    z,y,x = unzip(spots)
    res = list(zip(y,x))

    return res

def from_label_get_centeroidscoords(label):
    """Returns"""

    check_parameter(label = (np.ndarray))

    properties_dic = regionprops_table(label, properties= ["label","centroid"])
    Centroid = properties_dic
    return Centroid


