import bigfish.stack as stack
import numpy as np

def remove_mean_gaussian_background(image, sigma = 5):
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

    stack.check_parameter(image = (np.ndarray), sigma = (int, float))
    dim = image.ndim
    if dim == 2 :
        image_no_background = stack.remove_background_gaussian(image, sigma)
    
    elif dim == 3:
        mean_gaussian_background = stack.gaussian_filter(image, sigma)
        mean_gaussian_background = np.mean(mean_gaussian_background, axis= 0)
        image_no_background = np.subtract(image, mean_gaussian_background)
        image_no_background[image_no_background < 0] = 0
        image_no_background = np.round(image_no_background)
        image_no_background = np.array(image_no_background, dtype= image.dtype)
    
    else : raise Exception("Incorrect image dimension. Dimension should be either 2 or 3.")   

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

