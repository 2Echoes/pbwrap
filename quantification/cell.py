import numpy as np
from bigfish.stack import mean_projection, maximum_projection, check_parameter


from .utils import unzip


def mean_signal(cell, channel, projtype = 'mip') :
    """Returns mean signal from channel in computed cell. Mean signal is computed from 2D cell, so channel is projected z-wise according to projtype.
    
        Parameters
        ----------
            cell : dict
                Dictionary computed from bigFish.multistack.extract_cell
            channel : np.ndarray
                Channel from which intensity is computed
            projtype : str
                can either be mip or mean.

        Returns
        -------
            mean_sig : float
        
    """
    
    min_y, min_x, max_y, max_x = cell["bbox"]
    channel_cropped = channel[:, min_y:max_y, min_x:max_x]
    
    if projtype == 'mip' : 
        channel_crop_proj = maximum_projection(channel_cropped)
    elif projtype == 'mean' :
        channel_crop_proj = mean_projection(channel_cropped)
    min = channel_cropped.min()
    max = channel_crop_proj.max()
    mean = channel_crop_proj.mean()
    mean_sig = (mean-min)/(max-min)
    return mean_sig


def count_spots_in_mask(spots, mask) :
    """
    Parameters
    ----------

    """

    check_parameter(spots = (list), mask = np.ndarray)
    dim = len(spots[0])
    if dim != mask.ndim : raise ValueError("Number of coordinates whithin spots doesn't match mask's dimension.")
    if dim not in [2,3] : raise ValueError("Only 2D or 3D spots are supported.")
    spot_array = np.zeros_like(mask)
    
    if dim == 2 :
        line, col = unzip(spots)
        spot_array[line,col] = 1
    if dim == 3 : 
        plane,line,col = unzip(spots)
        spot_array[plane,line,col] = 1

    count_array = np.logical_and(spot_array, mask)
    print(count_array)
    count = count_array.sum()

    return count
