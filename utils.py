import numpy as np
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter
import time

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