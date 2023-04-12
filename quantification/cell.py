import numpy as np
from bigfish.stack import mean_projection, maximum_projection, check_parameter
from ..utils import rm_value_from_list
from .utils import unzip


def mean_signal(cell, channel, projtype = 'mip') :
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
    
    nucleus_mask = cell["nuc_mask"]
    metrics = compute_signalmetrics(channel_crop_proj, nucleus_mask)
    min = metrics["1 percentile"]
    max = metrics["99 percentile"]
    mean = metrics["mean"]
    mean_sig = (mean)/(max-min)
    return mean_sig


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
    #elif dim > 3 : raise Warning("spots are more than 3 dimensional, only first 3 axis will be computed.")

    
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
