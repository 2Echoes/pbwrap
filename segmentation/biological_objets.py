import numpy as np
import bigfish.stack as stack
import bigfish.segmentation as seg
import scipy.ndimage as ndi
from ..preprocessing import dog_filter, variance_filter
from .utils import get_histogramm_highest_varation_value, unstack_slices, from2Dlabel_to3Dlabel, auto_LoG_threshold
from .custom_functions import thresholding
from ..utils import compute_voxel_volume, gaussian_kernel_size
from ..errors import PbodySegmentationError
from bigfish.stack import check_array, check_parameter
from skimage.morphology import remove_small_objects


def pbody_segmentation(egfp, sigma = 2, threshold= None, thresh_penalty= 1, small_object_size= None, fill_holes= True) :
    """Performs Pbody segmentation on 2D or 3D egfp numpy array.
        Apply a log filter to the image which kernel is defined by sigma.
        Then a threshold is applied, if none is given compute automatic threshold from highest variation point in array histogram.
        if peaks_min_distance other than 0 is given, performs a watershed segmentation.
    
    Parameters
    ----------
        egfp : np.ndarray(y,x) or (z,y,x)
        beta : scalar
            Parameter representing the difficulty to cross region with high intensity gradient.
            
    Returns
    -------
        Pbody_label : np.ndarray(y,x) or (z,y,x)
    """

    check_parameter(egfp = (np.ndarray))
    check_array(egfp, ndim= [2,3])
    dim = egfp.ndim
    #Segmentation
    mask = stack.log_filter(egfp,sigma)
    if dim == 2 :
        if threshold == None : 
            threshold = get_histogramm_highest_varation_value(mask)
            threshold *= thresh_penalty
        mask = seg.thresholding(mask, threshold)
        mask = seg.clean_segmentation(mask, small_object_size=small_object_size, fill_holes=fill_holes)
    else : 
        for z in range(0,egfp.shape[0]):
            if threshold == None : 
                threshold = get_histogramm_highest_varation_value(mask[z])
                threshold *= thresh_penalty
            mask[z] = seg.thresholding(mask[z], threshold)
            mask[z] = seg.clean_segmentation(np.array(mask[z],dtype= bool), small_object_size=small_object_size, fill_holes=fill_holes)

    #labelling
    if dim == 2 : 
        egfp_label = seg.label_instances(mask)
        egfp_label = np.array(egfp_label, dtype= np.int64)

    else :
        slice_list = unstack_slices(mask)
        egfp_label = from2Dlabel_to3Dlabel(slice_list, maximal_distance= 10)



    if len(egfp_label) == 0 : raise PbodySegmentationError("No pbody was segmentated.")
    return egfp_label


def centrosome_segmentation_candidate_regions(image, centrosome_size, voxel_size, gaussian_fit = 'FWTM', threshold_penalty = 1, DoG_filter_method=True, LoG_filter_method=False, variance_method=True,) :
    """
    1st built for centrosome pipeline, centrosome segmentation.
    This function segments candidate regions for centrosme segmentation. More precise segmentation is then performed when label is processed cell by cell.

    Parameters
    ----------

    Returns
    -------
    labels : np.ndarray(z,y,x)
        3D labelled image, contains artifacts to be cleaned, and potentially more than 1 or 1 pair of centrosome per cell.
    """


    centrosome_voxel_ratio = [c_size / v_size for c_size, v_size in zip(centrosome_size, voxel_size)]
    centrosome_pixel_volume = compute_voxel_volume(centrosome_voxel_ratio)
    
    if DoG_filter_method and not LoG_filter_method :
        im_dog = dog_filter(image, kernel_size_1=1, kernel_size_2=3, voxel_size=voxel_size)
        threshold = np.percentile(im_dog, 99.96)
        mask = thresholding(im_dog, threshold)
        
    elif LoG_filter_method and not DoG_filter_method :
        im_log = stack.log_filter(image, sigma= centrosome_gaussian_size)
        centrosome_gaussian_size = gaussian_kernel_size(centrosome_size, voxel_size, width= gaussian_fit)
        flat_log = im_log.flatten()
        threshold = auto_LoG_threshold(flat_log) * threshold_penalty
        mask = thresholding(im_log, threshold)

    else : raise ValueError("Choose either DoG_filter_method (recommended for Alexa 647 fluorophore) or LoG_filter_method (recommended for EGFP)")
    
    if variance_method :
        im_var = variance_filter(image, (1,1,1))
        threshold = np.percentile(im_dog, 99.96)
        mask_ = thresholding(im_dog, threshold)

    
    mask = remove_small_objects(mask, round(centrosome_pixel_volume / 2**3)) # removing objects smaller than the half of expected centrosome
    labels,_ = ndi.label(mask)

    return labels