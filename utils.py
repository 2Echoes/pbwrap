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