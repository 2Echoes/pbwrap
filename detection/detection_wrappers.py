import bigfish.stack as stack
import bigfish.detection as detection
from bigfish.detection.spot_detection import local_maximum_detection, get_object_radius_pixel, _get_candidate_thresholds, spots_thresholding, _get_spot_counts
from ..errors import NoSpotError
import numpy as np
from types import GeneratorType


def spot_decomposition_nobckgrndrmv(image, spots, spot_radius, voxel_size_nm, alpha= 0.5, beta= 1):
    """ Basically same function as bigfish.detection.decompose_dense but without the remove background gaussian.
    
    Detect dense and bright regions with potential clustered spots and
    simulate a more realistic number of spots in these regions.

    #. We build a reference spot by aggregating predetected spots.
    #. We fit gaussian parameters on the reference spots.
    #. We detect dense regions to decompose.
    #. We simulate as many gaussians as possible in the candidate regions.

    Parameters
    ----------
    image : np.ndarray
        Image with shape (z, y, x) or (y, x).
    spots : np.ndarray
        Coordinate of the spots with shape (nb_spots, 3) or (nb_spots, 2)
        for 3-d or 2-d images respectively.
    voxel_size : int, float, Tuple(int, float) or List(int, float)
        Size of a voxel, in nanometer. One value per spatial dimension (zyx or
        yx dimensions). If it's a scalar, the same value is applied to every
        dimensions.
    spot_radius : int, float, Tuple(int, float) or List(int, float)
        Radius of the spot, in nanometer. One value per spatial dimension (zyx
        or yx dimensions). If it's a scalar, the same radius is applied to
        every dimensions.
    kernel_size : int, float, Tuple(float, int), List(float, int) or None
        Standard deviation used for the gaussian kernel (one for each
        dimension), in pixel. If it's a scalar, the same standard deviation is
        applied to every dimensions. If None, we estimate the kernel size from
        'spot_radius', 'voxel_size' and 'gamma'
    alpha : int or float
        Intensity percentile used to compute the reference spot, between 0
        and 1. The higher, the brighter are the spots simulated in the dense
        regions. Consequently, a high intensity score reduces the number of
        spots added. Default is 0.5, meaning the reference spot considered is
        the median spot.
    beta : int or float
        Multiplicative factor for the intensity threshold of a dense region.
        Default is 1. Threshold is computed with the formula:

        .. math::
            \\mbox{threshold} = \\beta * \\mbox{max(median spot)}

        With :math:`\\mbox{median spot}` the median value of all detected spot
        signals.
    gamma : int or float
        Multiplicative factor use to compute the gaussian kernel size:

        .. math::
            \\mbox{kernel size} = \\frac{\\gamma * \\mbox{spot radius}}{\\mbox{
            voxel size}}

        We perform a large gaussian filter with such scale to estimate image
        background and remove it from original image. A large gamma increases
        the scale of the gaussian filter and smooth the estimated background.
        To decompose very large bright areas, a larger gamma should be set.

    Notes
    -----
    If ``gamma = 0`` and ``kernel_size = None``, image is not denoised.

    Returns
    -------
    spots_postdecom : np.ndarray
        Coordinate of the spots detected, with shape (nb_spots, 3) or
        (nb_spots, 2). One coordinate per dimension (zyx or yx coordinates).

    """
    ndim = image.ndim


    #Spot decomposition    
    if len(spots) == 0 : return spots
    # reference spot
    reference_spot = detection.build_reference_spot(
        image=image,
        spots=spots,
        voxel_size=voxel_size_nm, 
        spot_radius= spot_radius,
        alpha=alpha)

    # fit a gaussian function on the reference spot
    sigma_z, sigma_yx, amplitude, background = detection.modelize_spot(
        reference_spot=reference_spot, 
        voxel_size= voxel_size_nm, 
        spot_radius= spot_radius)

    # detect dense regions
    regions_to_decompose, spots_out_regions, region_size = detection.get_dense_region(
        image= image, 
        spots= spots,
        voxel_size=voxel_size_nm,
        spot_radius= spot_radius,
        beta= beta)

    # precompute gaussian function values
    max_grid = max(200, region_size + 1)
    precomputed_gaussian = detection.precompute_erf(
        ndim= ndim,
        voxel_size=voxel_size_nm,
        sigma=(sigma_z, sigma_yx, sigma_yx),
        max_grid=max_grid)

    # simulate gaussian mixtures
    spots_in_regions, _ = detection.simulate_gaussian_mixture(
        image= image,
        candidate_regions=regions_to_decompose,
        voxel_size= voxel_size_nm,
        sigma=(sigma_z, sigma_yx, sigma_yx),
        amplitude=amplitude,
        background=background,
        precomputed_gaussian=precomputed_gaussian)

    spots_postdecomp = np.concatenate((spots_out_regions, spots_in_regions[:, :3]), axis=0)

    return spots_postdecomp




def mono_threshold_detect_spots(image,
        threshold=None,
        remove_duplicate=True,
        return_threshold=False,
        voxel_size=None,
        spot_radius=None,
        log_kernel_size=None,
        minimum_distance=None,
        threshold_factor= None):
    """
    Pbwrap
    ------
    Uses bigfish.detection function
    Can only take 1 image in input. 
    Add threshold factor which is multiplied to the auto_threshold from bigfish.
    For example to apply 20% penalty on automatic threshold threshold_factor should equal 1.20.
    
    Apply LoG filter followed by a Local Maximum algorithm to detect spots
    in a 2-d or 3-d image.

    #. We smooth the image with a LoG filter.
    #. We apply a multidimensional maximum filter.
    #. A pixel which has the same value in the original and filtered images
       is a local maximum.
    #. We remove local peaks under a threshold.
    #. We keep only one pixel coordinate per detected spot.

    Parameters
    ----------
    images : List[np.ndarray] or np.ndarray
        Image (or list of images) with shape (z, y, x) or (y, x). If several
        images are provided, the same threshold is applied.
    threshold : int, float or None
        A threshold to discriminate relevant spots from noisy blobs. If None,
        optimal threshold is selected automatically. If several images are
        provided, one optimal threshold is selected for all the images.
    remove_duplicate : bool
        Remove potential duplicate coordinates for the same spots. Slow the
        running.
    return_threshold : bool
        Return the threshold used to detect spots.
    voxel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of a voxel, in nanometer. One value per spatial dimension (zyx or
        yx dimensions). If it's a scalar, the same value is applied to every
        dimensions. Not used if 'log_kernel_size' and 'minimum_distance' are
        provided.
    spot_radius : int, float, Tuple(int, float), List(int, float) or None
        Radius of the spot, in nanometer. One value per spatial dimension (zyx
        or yx dimensions). If it's a scalar, the same radius is applied to
        every dimensions. Not used if 'log_kernel_size' and 'minimum_distance'
        are provided.
    log_kernel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of the LoG kernel. It equals the standard deviation (in pixels)
        used for the gaussian kernel (one for each dimension). One value per
        spatial dimension (zyx or yx dimensions). If it's a scalar, the same
        standard deviation is applied to every dimensions. If None, we estimate
        it with the voxel size and spot radius.
    minimum_distance : int, float, Tuple(int, float), List(int, float) or None
        Minimum distance (in pixels) between two spots we want to be able to
        detect separately. One value per spatial dimension (zyx or yx
        dimensions). If it's a scalar, the same distance is applied to every
        dimensions. If None, we estimate it with the voxel size and spot
        radius.

    Returns
    -------
    spots : List[np.ndarray] or np.ndarray, np.int64
        Coordinates (or list of coordinates) of the spots with shape
        (nb_spots, 3) or (nb_spots, 2), for 3-d or 2-d images respectively.
    threshold : int or float
        Threshold used to discriminate spots from noisy blobs.

    """
    stack.check_parameter(threshold_factor= (int, float))
    ndim = image.ndim
    # check consistency between parameters - detection with voxel size and
    # spot radius
    if (voxel_size is not None and spot_radius is not None
            and log_kernel_size is None and minimum_distance is None):
        if isinstance(voxel_size, (tuple, list)):
            if len(voxel_size) != ndim:
                raise ValueError("'voxel_size' must be a scalar or a sequence "
                                 "with {0} elements.".format(ndim))
        else:
            voxel_size = (voxel_size,) * ndim
        if isinstance(spot_radius, (tuple, list)):
            if len(spot_radius) != ndim:
                raise ValueError("'spot_radius' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            spot_radius = (spot_radius,) * ndim
        log_kernel_size = get_object_radius_pixel(
            voxel_size_nm=voxel_size,
            object_radius_nm=spot_radius,
            ndim=ndim)
        minimum_distance = get_object_radius_pixel(
            voxel_size_nm=voxel_size,
            object_radius_nm=spot_radius,
            ndim=ndim)

    # check consistency between parameters - detection with kernel size and
    # minimal distance
    elif (voxel_size is None and spot_radius is None
          and log_kernel_size is not None and minimum_distance is not None):
        if isinstance(log_kernel_size, (tuple, list)):
            if len(log_kernel_size) != ndim:
                raise ValueError("'log_kernel_size' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            log_kernel_size = (log_kernel_size,) * ndim
        if isinstance(minimum_distance, (tuple, list)):
            if len(minimum_distance) != ndim:
                raise ValueError("'minimum_distance' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            minimum_distance = (minimum_distance,) * ndim

    # check consistency between parameters - detection in priority with kernel
    # size and minimal distance
    elif (voxel_size is not None and spot_radius is not None
          and log_kernel_size is not None and minimum_distance is not None):
        if isinstance(log_kernel_size, (tuple, list)):
            if len(log_kernel_size) != ndim:
                raise ValueError("'log_kernel_size' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            log_kernel_size = (log_kernel_size,) * ndim
        if isinstance(minimum_distance, (tuple, list)):
            if len(minimum_distance) != ndim:
                raise ValueError("'minimum_distance' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            minimum_distance = (minimum_distance,) * ndim

    # missing parameters
    else:
        raise ValueError("One of the two pairs of parameters ('voxel_size', "
                         "'spot_radius') or ('log_kernel_size', "
                         "'minimum_distance') should be provided.")

    
    #If needed compute bigfish auto threshold
    if  threshold_factor != None:
        
        if threshold == None:
            image_filtered = stack.log_filter(image, log_kernel_size) #TODO should not stay None
            mask_local_max = local_maximum_detection(image_filtered, minimum_distance)
            threshold = detection.automated_threshold_setting(image_filtered, mask_local_max)
        
        threshold *= threshold_factor
    
    spots,threshold = detection.detect_spots(images=image, threshold=threshold, remove_duplicate= remove_duplicate, return_threshold= True, voxel_size=voxel_size, spot_radius=spot_radius, log_kernel_size=log_kernel_size, minimum_distance=minimum_distance)
    
    if len(spots) == 0 :
        raise NoSpotError("No spots were detected.")

    if return_threshold : return spots, threshold
    return spots



##################################################################################################################
"""
Below code was directly modified from the bigfish package : https://github.com/fish-quant/big-fish

BSD 3-Clause License
Copyright © 2020, Arthur Imbert
All rights reserved.
"""

def detect_spots(
        images,
        ndim,
        threshold=None,
        threshold_penalty: float = None,
        remove_duplicate=True,
        return_threshold=False,
        only_compute_threshold = False,
        voxel_size=None,
        spot_radius=None,
        log_kernel_size=None,
        minimum_distance=None,
        crop_zstack = None,
        show_advancement= False):

    """
    Pbwrap : In addition to original code we added : 
        threshold_penalty : float-like
            float which is multiplied with the automatic threshold setting.
        crop_zstack : list, tuple
         Crop images along z axis for Threshold calculation but NOT for spots detection.

    Apply LoG filter followed by a Local Maximum algorithm to detect spots
    in a 2-d or 3-d image.

    #. We smooth the image with a LoG filter.
    #. We apply a multidimensional maximum filter.
    #. A pixel which has the same value in the original and filtered images
       is a local maximum.
    #. We remove local peaks under a threshold.
    #. We keep only one pixel coordinate per detected spot.

    Parameters
    ----------
    images : List[np.ndarray] or np.ndarray
        Image (or list of images) with shape (z, y, x) or (y, x). If several
        images are provided, the same threshold is applied.
    threshold : int, float or None
        A threshold to discriminate relevant spots from noisy blobs. If None,
        optimal threshold is selected automatically. If several images are
        provided, one optimal threshold is selected for all the images.
    threshold_penalty : float
        When threshold is computed automatically this float is multiplied with 
        the computed threshold before application on spot detection.
    remove_duplicate : bool
        Remove potential duplicate coordinates for the same spots. Slow the
        running.
    return_threshold : bool
        Return the threshold used to detect spots.
    voxel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of a voxel, in nanometer. One value per spatial dimension (zyx or
        yx dimensions). If it's a scalar, the same value is applied to every
        dimensions. Not used if 'log_kernel_size' and 'minimum_distance' are
        provided.
    spot_radius : int, float, Tuple(int, float), List(int, float) or None
        Radius of the spot, in nanometer. One value per spatial dimension (zyx
        or yx dimensions). If it's a scalar, the same radius is applied to
        every dimensions. Not used if 'log_kernel_size' and 'minimum_distance'
        are provided.
    log_kernel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of the LoG kernel. It equals the standard deviation (in pixels)
        used for the gaussian kernel (one for each dimension). One value per
        spatial dimension (zyx or yx dimensions). If it's a scalar, the same
        standard deviation is applied to every dimensions. If None, we estimate
        it with the voxel size and spot radius.
    minimum_distance : int, float, Tuple(int, float), List(int, float) or None
        Minimum distance (in pixels) between two spots we want to be able to
        detect separately. One value per spatial dimension (zyx or yx
        dimensions). If it's a scalar, the same distance is applied to every
        dimensions. If None, we estimate it with the voxel size and spot
        radius.

    Returns
    -------
    spots : List[np.ndarray] or np.ndarray, np.int64
        Coordinates (or list of coordinates) of the spots with shape
        (nb_spots, 3) or (nb_spots, 2), for 3-d or 2-d images respectively.
    threshold : int or float
        Threshold used to discriminate spots from noisy blobs.

    """
    # check parameters
    stack.check_parameter(
        threshold=(int, float, type(None)),
        threshold_penalty= (int,float,type(None)),
        remove_duplicate=bool,
        return_threshold=bool,
        only_compute_threshold = bool,
        voxel_size=(int, float, tuple, list, type(None)),
        spot_radius=(int, float, tuple, list, type(None)),
        log_kernel_size=(int, float, tuple, list, type(None)),
        minimum_distance=(int, float, tuple, list, type(None)))

    # if one image is provided we enlist it
    if type(crop_zstack) != type(None):
        if not isinstance(crop_zstack, (list,tuple)) : raise TypeError("Wrong type for crop_ztack it should be None, tuple or list. Is is {0}".format(type(crop_zstack)))
        if len(crop_zstack) != 2 : raise ValueError("crop_zstack is expected to have 2 elements. It has {0}".format(len(crop_zstack)))
        if not (isinstance(crop_zstack[0],int) and isinstance(crop_zstack[1], int)) : raise TypeError("crop_zstack elements should be ints.")
        if crop_zstack[1] <= crop_zstack[0] : raise ValueError("crop zstack should be [first crop, last_crop[ last_crop can't be smaller than or equals first crop.")

    if not (isinstance(images, list) or isinstance(images, GeneratorType)):
        stack.check_array(
            images,
            ndim=[2, 3],
            dtype=[np.uint8, np.uint16, np.float32, np.float64])
        ndim = images.ndim
        images = [images]
        is_list = False
    elif isinstance(images, list):
        ndim = None
        for i, image in enumerate(images):
            stack.check_array(
                image,
                ndim=[2, 3],
                dtype=[np.uint8, np.uint16, np.float32, np.float64])
            if i == 0:
                ndim = image.ndim
            else:
                if ndim != image.ndim:
                    raise ValueError("Provided images should have the same "
                                     "number of dimensions.")
        is_list = True
    else : is_list = True

    # check consistency between parameters - detection with voxel size and
    #return threshold
    if only_compute_threshold == True and return_threshold == False :
        raise ValueError("if return_threshold is False, only_compute_threshold should not be set to True.")
    # spot radius
    if (voxel_size is not None and spot_radius is not None
            and log_kernel_size is None and minimum_distance is None):
        if isinstance(voxel_size, (tuple, list)):
            if len(voxel_size) != ndim:
                raise ValueError("'voxel_size' must be a scalar or a sequence "
                                 "with {0} elements.".format(ndim))
        else:
            voxel_size = (voxel_size,) * ndim
        if isinstance(spot_radius, (tuple, list)):
            if len(spot_radius) != ndim:
                raise ValueError("'spot_radius' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            spot_radius = (spot_radius,) * ndim
        log_kernel_size = get_object_radius_pixel(
            voxel_size_nm=voxel_size,
            object_radius_nm=spot_radius,
            ndim=ndim)
        minimum_distance = get_object_radius_pixel(
            voxel_size_nm=voxel_size,
            object_radius_nm=spot_radius,
            ndim=ndim)

    # check consistency between parameters - detection with kernel size and
    # minimal distance
    elif (voxel_size is None and spot_radius is None
          and log_kernel_size is not None and minimum_distance is not None):
        if isinstance(log_kernel_size, (tuple, list)):
            if len(log_kernel_size) != ndim:
                raise ValueError("'log_kernel_size' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            log_kernel_size = (log_kernel_size,) * ndim
        if isinstance(minimum_distance, (tuple, list)):
            if len(minimum_distance) != ndim:
                raise ValueError("'minimum_distance' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            minimum_distance = (minimum_distance,) * ndim

    # check consistency between parameters - detection in priority with kernel
    # size and minimal distance
    elif (voxel_size is not None and spot_radius is not None
          and log_kernel_size is not None and minimum_distance is not None):
        if isinstance(log_kernel_size, (tuple, list)):
            if len(log_kernel_size) != ndim:
                raise ValueError("'log_kernel_size' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            log_kernel_size = (log_kernel_size,) * ndim
        if isinstance(minimum_distance, (tuple, list)):
            if len(minimum_distance) != ndim:
                raise ValueError("'minimum_distance' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            minimum_distance = (minimum_distance,) * ndim

    # missing parameters
    else:
        raise ValueError("One of the two pairs of parameters ('voxel_size', "
                         "'spot_radius') or ('log_kernel_size', "
                         "'minimum_distance') should be provided.")

    # detect spots
    if return_threshold:
        if only_compute_threshold :
            threshold = _detect_spots_from_images(
                images,
                threshold=threshold,
                threshold_penalty= threshold_penalty,
                remove_duplicate=remove_duplicate,
                return_threshold=return_threshold,
                only_compute_threshold = only_compute_threshold,
                log_kernel_size=log_kernel_size,
                min_distance=minimum_distance,
                crop_zstack=crop_zstack,
                show_advancement=show_advancement)
        else :
            spots, threshold = _detect_spots_from_images(
                images,
                threshold=threshold,
                threshold_penalty= threshold_penalty,
                remove_duplicate=remove_duplicate,
                return_threshold=return_threshold,
                only_compute_threshold = only_compute_threshold,
                log_kernel_size=log_kernel_size,
                min_distance=minimum_distance,
                crop_zstack=crop_zstack,
                show_advancement=show_advancement)
    else:
        spots = _detect_spots_from_images(
            images,
            threshold=threshold,
            threshold_penalty= threshold_penalty,
            remove_duplicate=remove_duplicate,
            return_threshold=return_threshold,
            log_kernel_size=log_kernel_size,
            min_distance=minimum_distance,
            crop_zstack=crop_zstack,
            show_advancement=show_advancement)

    # format results
    if only_compute_threshold : return threshold

    if not is_list:
        spots = spots[0]

    # return threshold or not
    if return_threshold:
        return spots, threshold
    else:
        return spots





def _detect_spots_from_images(
        images,
        threshold=None,
        threshold_penalty:float = None,
        remove_duplicate=True,
        return_threshold=False,
        only_compute_threshold= False,
        log_kernel_size=None,
        min_distance=None,
        crop_zstack= None,
        show_advancement=False):
    """Apply LoG filter followed by a Local Maximum algorithm to detect spots
    in a 2-d or 3-d image.

    #. We smooth the image with a LoG filter.
    #. We apply a multidimensional maximum filter.
    #. A pixel which has the same value in the original and filtered images
       is a local maximum.
    #. We remove local peaks under a threshold.
    #. We keep only one pixel coordinate per detected spot.

    Parameters
    ----------
    images : List[np.ndarray]
        List of images with shape (z, y, x) or (y, x). The same threshold is
        applied to every images.
    threshold : float or int
        A threshold to discriminate relevant spots from noisy blobs. If None,
        optimal threshold is selected automatically. If several images are
        provided, one optimal threshold is selected for all the images.
    remove_duplicate : bool
        Remove potential duplicate coordinates for the same spots. Slow the
        running.
    return_threshold : bool
        Return the threshold used to detect spots.
    log_kernel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of the LoG kernel. It equals the standard deviation (in pixels)
        used for the gaussian kernel (one for each dimension). One value per
        spatial dimension (zyx or yx dimensions). If it's a scalar, the same
        standard deviation is applied to every dimensions. If None, we estimate
        it with the voxel size and spot radius.
    min_distance : int, float, Tuple(int, float), List(int, float) or None
        Minimum distance (in pixels) between two spots we want to be able to
        detect separately. One value per spatial dimension (zyx or yx
        dimensions). If it's a scalar, the same distance is applied to every
        dimensions. If None, we estimate it with the voxel size and spot
        radius.

    Returns
    -------
    all_spots : List[np.ndarray], np.int64
        List of spot coordinates with shape (nb_spots, 3) or (nb_spots, 2),
        for 3-d or 2-d images respectively.
    threshold : int or float
        Threshold used to discriminate spots from noisy blobs.

    """
    if threshold_penalty == None : threshold_penalty = 1
    # initialization
    n = 0

    # apply LoG filter and find local maximum
    images_filtered = []
    pixel_values = []
    masks = []
    for image in images:
        n += 1
        # filter image
        image_filtered = stack.log_filter(image, log_kernel_size)
        images_filtered.append(image_filtered)

        # get pixels value
        if crop_zstack != None :
            pixel_values += list(image_filtered[crop_zstack[0]:crop_zstack[1]].ravel())
        else : 
            pixel_values += list(image_filtered.ravel())

        # find local maximum
        mask_local_max = local_maximum_detection(image_filtered, min_distance)
        masks.append(mask_local_max)

    # get optimal threshold if necessary based on all the images
    if threshold is None:

        # get threshold values we want to test
        thresholds = _get_candidate_thresholds(pixel_values)

        # get spots count and its logarithm
        all_value_spots = []
        minimum_threshold = float(thresholds[0])
        for i in range(n):
            image_filtered = images_filtered[i]
            mask_local_max = masks[i]
            spots, mask_spots = spots_thresholding(
                image_filtered, mask_local_max,
                threshold=minimum_threshold,
                remove_duplicate=False)
            value_spots = image_filtered[mask_spots]
            all_value_spots.append(value_spots)
        all_value_spots = np.concatenate(all_value_spots)
        thresholds, count_spots = _get_spot_counts(thresholds, all_value_spots)

        # select threshold where the kink of the distribution is located
        if count_spots.size > 0:
            threshold, _, _ = get_breaking_point(thresholds, count_spots)
            threshold *= threshold_penalty
            if only_compute_threshold : return threshold

    # detect spots
    all_spots = []
    for i in range(n):
        if show_advancement : print("image {0}/{1}".format(i+1,n))

        # get images and masks
        image_filtered = images_filtered[i]
        mask_local_max = masks[i]

        # detection
        spots, _ = spots_thresholding(
            image_filtered, mask_local_max, threshold, remove_duplicate)
        all_spots.append(spots)

    # return threshold or not
    if return_threshold:
        return all_spots, threshold
    else:
        return all_spots
    




def iter_detect_spots(
        images,
        ndim: int,
        threshold=None,
        threshold_penalty: float = None,
        remove_duplicate=True,
        return_threshold=False,
        voxel_size=None,
        spot_radius=None,
        log_kernel_size=None,
        minimum_distance=None):

    """
    Pbwrap : In addition to original code we added a parameter : threshold_penalty : float which is multiplied with the automatic threshold setting.

    Apply LoG filter followed by a Local Maximum algorithm to detect spots
    in a 2-d or 3-d image.

    #. We smooth the image with a LoG filter.
    #. We apply a multidimensional maximum filter.
    #. A pixel which has the same value in the original and filtered images
       is a local maximum.
    #. We remove local peaks under a threshold.
    #. We keep only one pixel coordinate per detected spot.

    Parameters
    ----------
    images : List[np.ndarray] or np.ndarray
        Image (or list of images) with shape (z, y, x) or (y, x). If several
        images are provided, the same threshold is applied.
    threshold : int, float or None
        A threshold to discriminate relevant spots from noisy blobs. If None,
        optimal threshold is selected automatically. If several images are
        provided, one optimal threshold is selected for all the images.
    threshold_penalty : float
        When threshold is computed automatically this float is multiplied with 
        the computed threshold before application on spot detection.
    remove_duplicate : bool
        Remove potential duplicate coordinates for the same spots. Slow the
        running.
    return_threshold : bool
        Return the threshold used to detect spots.
    voxel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of a voxel, in nanometer. One value per spatial dimension (zyx or
        yx dimensions). If it's a scalar, the same value is applied to every
        dimensions. Not used if 'log_kernel_size' and 'minimum_distance' are
        provided.
    spot_radius : int, float, Tuple(int, float), List(int, float) or None
        Radius of the spot, in nanometer. One value per spatial dimension (zyx
        or yx dimensions). If it's a scalar, the same radius is applied to
        every dimensions. Not used if 'log_kernel_size' and 'minimum_distance'
        are provided.
    log_kernel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of the LoG kernel. It equals the standard deviation (in pixels)
        used for the gaussian kernel (one for each dimension). One value per
        spatial dimension (zyx or yx dimensions). If it's a scalar, the same
        standard deviation is applied to every dimensions. If None, we estimate
        it with the voxel size and spot radius.
    minimum_distance : int, float, Tuple(int, float), List(int, float) or None
        Minimum distance (in pixels) between two spots we want to be able to
        detect separately. One value per spatial dimension (zyx or yx
        dimensions). If it's a scalar, the same distance is applied to every
        dimensions. If None, we estimate it with the voxel size and spot
        radius.

    Returns
    -------
    spots : List[np.ndarray] or np.ndarray, np.int64
        Coordinates (or list of coordinates) of the spots with shape
        (nb_spots, 3) or (nb_spots, 2), for 3-d or 2-d images respectively.
    threshold : int or float
        Threshold used to discriminate spots from noisy blobs.

    """
    # check parameters
    stack.check_parameter(
        threshold=(int, float, type(None)),
        threshold_penalty= (int,float,type(None)),
        remove_duplicate=bool,
        return_threshold=bool,
        voxel_size=(int, float, tuple, list, type(None)),
        spot_radius=(int, float, tuple, list, type(None)),
        log_kernel_size=(int, float, tuple, list, type(None)),
        minimum_distance=(int, float, tuple, list, type(None)))
    
    is_list = True

    # check consistency between parameters - detection with voxel size and
    # spot radius
    if (voxel_size is not None and spot_radius is not None
            and log_kernel_size is None and minimum_distance is None):
        if isinstance(voxel_size, (tuple, list)):
            if len(voxel_size) != ndim:
                raise ValueError("'voxel_size' must be a scalar or a sequence "
                                 "with {0} elements.".format(ndim))
        else:
            voxel_size = (voxel_size,) * ndim
        if isinstance(spot_radius, (tuple, list)):
            if len(spot_radius) != ndim:
                raise ValueError("'spot_radius' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            spot_radius = (spot_radius,) * ndim
        log_kernel_size = get_object_radius_pixel(
            voxel_size_nm=voxel_size,
            object_radius_nm=spot_radius,
            ndim=ndim)
        minimum_distance = get_object_radius_pixel(
            voxel_size_nm=voxel_size,
            object_radius_nm=spot_radius,
            ndim=ndim)

    # check consistency between parameters - detection with kernel size and
    # minimal distance
    elif (voxel_size is None and spot_radius is None
          and log_kernel_size is not None and minimum_distance is not None):
        if isinstance(log_kernel_size, (tuple, list)):
            if len(log_kernel_size) != ndim:
                raise ValueError("'log_kernel_size' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            log_kernel_size = (log_kernel_size,) * ndim
        if isinstance(minimum_distance, (tuple, list)):
            if len(minimum_distance) != ndim:
                raise ValueError("'minimum_distance' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            minimum_distance = (minimum_distance,) * ndim

    # check consistency between parameters - detection in priority with kernel
    # size and minimal distance
    elif (voxel_size is not None and spot_radius is not None
          and log_kernel_size is not None and minimum_distance is not None):
        if isinstance(log_kernel_size, (tuple, list)):
            if len(log_kernel_size) != ndim:
                raise ValueError("'log_kernel_size' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            log_kernel_size = (log_kernel_size,) * ndim
        if isinstance(minimum_distance, (tuple, list)):
            if len(minimum_distance) != ndim:
                raise ValueError("'minimum_distance' must be a scalar or a "
                                 "sequence with {0} elements.".format(ndim))
        else:
            minimum_distance = (minimum_distance,) * ndim

    # missing parameters
    else:
        raise ValueError("One of the two pairs of parameters ('voxel_size', "
                         "'spot_radius') or ('log_kernel_size', "
                         "'minimum_distance') should be provided.")

    # detect spots
    if return_threshold:
        spots, threshold = iter_detect_spots_from_images(
            images,
            threshold=threshold,
            threshold_penalty= threshold_penalty,
            remove_duplicate=remove_duplicate,
            return_threshold=return_threshold,
            log_kernel_size=log_kernel_size,
            min_distance=minimum_distance)
    else:
        spots = iter_detect_spots_from_images(
            images,
            threshold=threshold,
            threshold_penalty= threshold_penalty,
            remove_duplicate=remove_duplicate,
            return_threshold=return_threshold,
            log_kernel_size=log_kernel_size,
            min_distance=minimum_distance)

    # return threshold or not
    if return_threshold:
        return spots, threshold
    else:
        return spots





def iter_detect_spots_from_images(
        images,
        threshold=None,
        threshold_penalty:float = None,
        remove_duplicate=True,
        return_threshold=False,
        log_kernel_size=None,
        min_distance=None):
    """Apply LoG filter followed by a Local Maximum algorithm to detect spots
    in a 2-d or 3-d image.

    #. We smooth the image with a LoG filter.
    #. We apply a multidimensional maximum filter.
    #. A pixel which has the same value in the original and filtered images
       is a local maximum.
    #. We remove local peaks under a threshold.
    #. We keep only one pixel coordinate per detected spot.

    Parameters
    ----------
    images : List[np.ndarray]
        List of images with shape (z, y, x) or (y, x). The same threshold is
        applied to every images.
    threshold : float or int
        A threshold to discriminate relevant spots from noisy blobs. If None,
        optimal threshold is selected automatically. If several images are
        provided, one optimal threshold is selected for all the images.
    remove_duplicate : bool
        Remove potential duplicate coordinates for the same spots. Slow the
        running.
    return_threshold : bool
        Return the threshold used to detect spots.
    log_kernel_size : int, float, Tuple(int, float), List(int, float) or None
        Size of the LoG kernel. It equals the standard deviation (in pixels)
        used for the gaussian kernel (one for each dimension). One value per
        spatial dimension (zyx or yx dimensions). If it's a scalar, the same
        standard deviation is applied to every dimensions. If None, we estimate
        it with the voxel size and spot radius.
    min_distance : int, float, Tuple(int, float), List(int, float) or None
        Minimum distance (in pixels) between two spots we want to be able to
        detect separately. One value per spatial dimension (zyx or yx
        dimensions). If it's a scalar, the same distance is applied to every
        dimensions. If None, we estimate it with the voxel size and spot
        radius.

    Returns
    -------
    all_spots : List[np.ndarray], np.int64
        List of spot coordinates with shape (nb_spots, 3) or (nb_spots, 2),
        for 3-d or 2-d images respectively.
    threshold : int or float
        Threshold used to discriminate spots from noisy blobs.

    """
    if threshold_penalty == None : threshold_penalty = 1
    # initialization
    n = 0

    # apply LoG filter and find local maximum
    pixel_values = []
    images_filtered = [] 
    for image in images : #Loading data into RAM
        image_filtered = stack.log_filter(image,log_kernel_size)
        images_filtered += [image_filtered]
        pixel_values += list(image_filtered.ravel())
    n = len(images_filtered)
    masks = (local_maximum_detection(image_filtered, min_distance) for image_filtered in images_filtered)
    #Re-generating images_filtered after masks computation

    # get optimal threshold if necessary based on all the images
    if threshold is None:

        # get threshold values we want to test
        thresholds = _get_candidate_thresholds(pixel_values)

        # get spots count and its logarithm
        minimum_threshold = float(thresholds[0])
        threshold_list = []
        for image_filtered, mask_local_max in zip(images_filtered, masks) :
            all_value_spots = []
            spots, mask_spots = spots_thresholding(
                image_filtered, mask_local_max,
                threshold=minimum_threshold,
                remove_duplicate=False)
            value_spots = image_filtered[mask_spots]
            all_value_spots.append(value_spots)
            all_value_spots = np.concatenate(all_value_spots)
            thresholds, count_spots = _get_spot_counts(thresholds, all_value_spots)

        # select threshold where the kink of the distribution is located
            if count_spots.size > 0:
                threshold, _, _ = get_breaking_point(thresholds, count_spots)
            threshold_list += [threshold]



        for i in range(n):
            all_value_spots = []
            image_filtered = images_filtered[i]
            mask_local_max = masks[i]
            spots, mask_spots = spots_thresholding(
                image_filtered, mask_local_max,
                threshold=minimum_threshold,
                remove_duplicate=False)
            value_spots = image_filtered[mask_spots]
            all_value_spots.append(value_spots)
            all_value_spots = np.concatenate(all_value_spots)
            thresholds, count_spots = _get_spot_counts(thresholds, all_value_spots)

        # select threshold where the kink of the distribution is located
            if count_spots.size > 0:
                threshold, _, _ = get_breaking_point(thresholds, count_spots)
            threshold_list += [threshold]
        threshold = np.median(threshold_list)
        threshold_list += [threshold]
        threshold = np.median(threshold_list)
        threshold *= threshold_penalty

    # detect spots
    all_spots = []
    for i in range(n):

        # get images and masks
        image_filtered = images_filtered[i]
        mask_local_max = masks[i]

        # detection
        spots, _ = spots_thresholding(
            image_filtered, mask_local_max, threshold, remove_duplicate)
        all_spots.append(spots)

    # return threshold or not
    if return_threshold:
        return all_spots, threshold
    else:
        return all_spots



def _compute_local_max_masks(images_filtered, min_distance) :
    """Return iterator of local maximum masks computed from log filtered images"""

    pixel_values = []
    for image_filtered in images_filtered :
        
        # get pixels value
        pixel_values += list(image_filtered.ravel())

        # find local maximum
        yield local_maximum_detection(image_filtered, min_distance)





def get_breaking_point(x, y):
    """Select the x-axis value where a L-curve has a kink.

    Assuming a L-curve from A to B, the 'breaking_point' is the more distant
    point to the segment [A, B].

    Parameters
    ----------
    x : np.array
        X-axis values.
    y : np.array
        Y-axis values.

    Returns
    -------
    breaking_point : float
        X-axis value at the kink location.
    x : np.array
        X-axis values.
    y : np.array
        Y-axis values.

    """
    # check parameters
    stack.check_array(
        x,
        ndim=1,
        dtype=[np.float32, np.float64, np.int32, np.int64])
    stack.check_array(
        y,
        ndim=1,
        dtype=[np.float32, np.float64, np.int32, np.int64])

    # select threshold where curve break
    slope = (y[-1] - y[0]) / len(y)
    y_grad = np.gradient(y)
    m = list(y_grad >= slope)
    j = m.index(False)
    m = m[j:]
    x = x[j:]
    y = y[j:]
    if True in m:
        i = m.index(True)
    else:
        i = -1
    breaking_point = float(x[i])

    return breaking_point, x, y