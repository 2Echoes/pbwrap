import numpy as np
from bigfish.stack import check_parameter
from .errors import DetectionTimeOutError

def detectiontimeout_handler(signum, frame):
    raise DetectionTimeOutError('Acquisition processing timeout.')


def _is_sameshape(*arrays) :
    """Returns True if numpy *arrays have all the same shape or if al lists have the same length.
    
    Parameters
    ----------
        *arrays : np.ndarray or list
            arrays on which shape check is performed.
    
    Returns
    -------
        res : bool
            True if all arrays have the same shape else False
    """
    ref_type = type(arrays[0])
    truth_table = []

    if ref_type == list :
        ref_length = len(arrays[0])
        for array in arrays :
            truth_table += [ref_length == len(array)]

    if ref_type == np.ndarray :
        ref_shape = arrays[0].shape
        for array in arrays:
            truth_table += [array.shape == ref_shape]
    
    res = all(truth_table)
    return res


def check_sameshape(*arrays) :
    """Raises exception if all *arrays have not the same shapeo r if al lists have not the same length.
    
    Parameters
    ----------
        *arrays : np.ndarray or list
            arrays on which shape check is performed.
    
    """

    for array in arrays :
        check_parameter(array = (np.ndarray, list))

    if not _is_sameshape(*arrays): raise Exception("Arrays do not have the same shape.")