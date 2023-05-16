import numpy as np
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter

def from_label_get_centeroidscoords(label):
    """Returns dict{"label", "centroid"}"""

    check_parameter(label = (np.ndarray))
    dim = label.ndim
    properties_dic = regionprops_table(label, properties= ["label","centroid"])
    centroid = properties_dic
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
    #obsolete use list.remove(value) instead
    index = get_elmtindex(value,lis)
    count = 0
    for i in index :
        if i-count == len(lis)-1 :
            lis = lis[:i-count]
            count += 1
        else:
            lis = lis[: i - count] + lis[i + 1 - count :]
            count += 1 
    return lis

def point_is_in_mask(point_coords: tuple, mask: np.ndarray):
    if len(point_coords) != mask.ndim : raise ValueError("Number of dimension betwenn point and mask doesn't match")
    if not isinstance(point_coords, (list,tuple)) : raise TypeError("wrong type for point_coords parameter")
    if not isinstance(mask, np.ndarray) : raise TypeError("wrong type for mask parameter")

    if mask.dtype != bool : mask = np.array(mask, dtype= bool)

    return mask[point_coords[0], point_coords[1]] == 1