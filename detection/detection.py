import bigfish.stack as stack
import bigfish.detection as detection
from bigfish.detection.spot_detection import local_maximum_detection, get_object_radius_pixel
import numpy as np

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




def detect_spots(image,
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
    stack.check_parameter(image = (np.ndarray), threshold_factor= (int, float))
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
    
    if return_threshold : return spots, threshold
    return spots