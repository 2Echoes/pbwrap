import numpy as np
import bigfish.stack as stack
import scipy.ndimage as ndi
from pbwrap.utils import compute_anisotropy_coef
from skimage.segmentation import random_walker, watershed
from skimage.feature import peak_local_max

def gaussian_threshold_segmentation(image : np.ndarray, sigma, percentile, voxel_size)-> 'np.ndarray[bool]' :
    """
    Naive segmentation : first gaussian blur of kernel sigma is applied, then a threshold keeping pixels higer than the image-percentile will be used.

    Parameters
    ----------
        image : np.ndarray
            Either 3D (z; y; x) or 2D (y; x) image.
        sigma : float-like
            Gaussian blur kernel size will be adapted for anisotropic image using voxel_size
        percentile : float
            Between 0+ and 100. Threshold keeps the higer 'percentile' pixels. --> 1 keeps no pixel, 25 keeps top 75% pixels.
        voxel_size : tuple
            same dim as image. represent scale factor between pixels and real distance.
    """

    stack.check_parameter(image = np.ndarray, sigma = (float, int), percentile = (float, int), voxel_size = (tuple, list))
    if len(voxel_size) != image.ndim : raise ValueError("Voxel size length must be equal to image dimension.")
    if percentile <= 0 or percentile > 100 : raise ValueError("percentile must be in ]0;1].")

    anisotropy_coef = compute_anisotropy_coef(voxel_size)
    sigma = [coef * sigma for coef in anisotropy_coef]

    image = stack.gaussian_filter(image, sigma=sigma)
    percentile_pixel = np.percentile(image, percentile)
    return image > percentile_pixel

def thresholding(image, threshold) :
    """
    Return thresholded image as a boolean array. threshold is strictly applied.
    """
    if not isinstance(threshold, (float,int)) : raise TypeError("threshold should be int or float, it is {0}".format(type(threshold)))
    if not isinstance(image, np.ndarray) :
        try :
            im = np.array(image, dtype= float)
        except Exception :
            raise ValueError("image parameter is not numpy arrray and couldn't be converted to float array.")
    else : im = image.copy()

    return im > threshold


def random_walker_segmentation(image, percentile_down = 99.5, percentile_up = 99.8, beta= 1):
    #TODO Unused
    """Performs random walker segmentation using scipy algorithm. The segmentation is performed by assigning seeds to element in the image.
    In our algorithm we're trying to seperate background from one type of object (mainly pbodies). We assign the seed 2 to pixels we know are in the background and seed 1 to pixels we know are in p-bodies.
    Pixel left at 0 are the pixels the random walker segment into group 1 (pbodies) or group 2 background.
    Afterwards, background is set back to 0 so output is a classical binary mask.

    Percentiles paremeters should be fine tuned to your image, default settings correspond to p-body seg using egfp channel.
    
    Parameters
    ----------
        image : np.ndarray
            3D or 2D image to segment.
        percentile_down : scalar (from 0 to 100)
            Percentile of pixel set into background. (group 2)
        percentile_up : scalar (from 0 to 100 and > percentile_down)
            100 - x highest percentile of pixels set into p-bodies
        beta : scalar
            Defines how hard it is to break intensity gradient during segmentation.

    Returns
    -------
        mask : np.ndarray(bool)

    """
    stack.check_parameter(image = (np.ndarray), percentile_down = (int, float), percentile_up= (int, float), beta= (int, float))
    if percentile_down < 0 or percentile_down > 100 : raise Exception("Percentile_down parameter should be in range 0-100.")
    if percentile_up < 0 or percentile_down > 100 : raise Exception("Percentile_up parameter should be in range 0-100.")
    if percentile_up < percentile_down : raise Exception("Percentile_up parameter should be larger than percentile_down to avoid conflit when attributing seeds.")

    seed = np.zeros_like(image)
    seed[image > np.percentile(image, percentile_up)] = 1
    seed[image < np.percentile(image, percentile_down)] = 2
    mask = random_walker(image, seed, beta=beta)
    mask[mask != 1] = 0
    mask = mask.astype(bool)
    
    return mask




def watershed_segmentation(image, label=None, peaks_min_distance= 3, inv_image = False ):
    #TODO : Add sampling [3,1,1] or [anisotropy, 1, 1] in case of 3D segmentation.
    #TODO : Unused
    """Performs watershed segmentation using scipy algorithm. 
    In the usual regions flooding thought process this algorithm uses local maxima as sources for the flooding. 
    
    Parameters
    ----------
        image : np.ndarray
            3D or 2D image. For optimal performance input image should be either a boolean image or a labelled image, for a grayscale image make sure to be restrictive enough on local maxima computation.
        peaks_min_distance : int
            Minimal distance (in  pixel) separating two maximum intensity peaks --> if d = 1 the maximum number of peaks is computed.
            
    Returns
    -------
        label : np.ndarray
            labelled image.
    """

    stack.check_parameter(image = (np.ndarray), peaks_min_distance = (int))
    if inv_image :
        image = np.invert(image)

    distance = ndi.distance_transform_edt(image)
    #distance = distance_transform(image, label)
    coords = peak_local_max(distance, footprint=np.ones((3, 3)), min_distance = peaks_min_distance) #labels = label
    mask_water = np.zeros(distance.shape, dtype=bool)
    mask_water[tuple(coords.T)] = True
    markers, _ = ndi.label(mask_water)
    res = watershed(-distance, markers, compactness= 100) #mask = label

    return res