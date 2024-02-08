import bigfish.stack as stack
import numpy as np

from pbwrap.utils import compute_anisotropy_coef

def remove_mean_gaussian_background(image, sigma = 5, voxel_size = (1,1,1)):
    """Removes background of the image using gaussian filter. If the input image is 3D the mean gaussian background is computed from all slices and then substract to the stack.
    
    Parameters
    ----------
        image : np.ndarray
    
        sigma : scalar
            Defines the gaussian filter intensity.
    Returns
    -------
        image_no_background : np.ndarray
            Image with background substracted.
    """

    if not len(voxel_size) == image.ndim : raise ValueError("Inconsistency between voxel_size length and image dimension.")
    anisotropy_coef = compute_anisotropy_coef(voxel_size=voxel_size)
    corrected_sigma = [sigma / anisotropy_i for anisotropy_i in anisotropy_coef]

    image_no_background = stack.remove_background_gaussian(image, corrected_sigma)

    return image_no_background



def get_first_infocus(image, score_threshold = 9):
    """Return index of first in focus slice from a 3D image
    
    Parameters
    ----------
        image : np.ndarray(z,y,x)
        
    Returns
    -------
        res : int
            returns index on z axis of the first in focus slice. returns -1 if score threshold is never reached
    """

    z = -1
    score = 0
    while score < score_threshold and z+1 < image.shape[0] :
        z +=1
        score = stack.compute_focus(image[z]).max()
    if z >= image.shape[0] : 
        raise Warning("No slice scored above the threshold (focus score)")
    else : return z

