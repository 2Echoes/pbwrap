import bigfish.stack as stack
import numpy as np
import scipy.ndimage as ndi
from pbwrap.utils import compute_anisotropy_coef

def _disk_unit_kernel(shape) :
    kernel = np.ones(shape= [2*axis -1 for axis in shape], dtype= int) #always has a center pixel because 2*axis-1 always an odd int
    center_coordinates = tuple([axis -1 for axis in shape])
    kernel[center_coordinates] = 0
    sampling = [1/(axis-1) if axis !=1 else 1 for axis in shape]
    distance_map = ndi.distance_transform_edt(kernel, sampling= sampling) #distance map to the center of the kernel, sampling is used for anisotropic kernels
    kernel[distance_map <= 1] = 1
    kernel[distance_map > 1] = 0

    return kernel

def _mean_kernel(shape) :
    mean_kernel = _disk_unit_kernel(shape).astype(float)
    mean_kernel /= mean_kernel.sum()

    return mean_kernel

def _unsigned_int_substract(arr1, arr2) :
    """
    if arr2[coord] > arr1[coord] then arr1[coord] - arr2[coord] is set to 0 instead of (2^n -1) - (arr1[coord] - arr2[coord])
    """

    
    if arr1.dtype != arr2.dtype : raise ValueError("arr1 and arr2 must have the same unsigned dtype.")
    else : dtype = arr1.dtype

    arr1 = arr1.astype(int)
    arr2 = arr2.astype(int)

    res = arr1 - arr2
    res[res < 0] = 0
    res.astype(dtype)
    return res


def mean_filter(image, shape) :
    """
    Only disk shaped kernel are supported.
    """
    if len(shape) != image.ndim : raise ValueError("Shape length and image dimension must match.")
    mean_kernel = _mean_kernel(shape)
    image = ndi.convolve(image, mean_kernel, mode= 'nearest')
    return image

def variance_filter(image, kernel_size) :
    stack.check_parameter(image = np.ndarray)
    squared_mean_image = np.power(
        mean_filter(image, shape=kernel_size)
        , 2)
    mean_squared_image = mean_filter(
        np.power(image, 2),
        shape=kernel_size)

    res = _unsigned_int_substract(squared_mean_image, mean_squared_image)

    return res

def dog_filter(image, kernel_size_1, kernel_size_2,  voxel_size) :
    """
    3D/2D Difference Of Gaussian (DOG) filter. Performs 2 3D gaussians with kernel_size 1 & 2 and performs the substraction : gauss1 - gauss2.
    """

    if not image.ndim == len(voxel_size) : raise ValueError("image dimension and voxel_size length must match.")

    anisotropy_coef = compute_anisotropy_coef(voxel_size=voxel_size)
    corrected_sigma_1 = [kernel_size_1 / anisotropy_i for anisotropy_i in anisotropy_coef]
    corrected_sigma_2 = [kernel_size_2 / anisotropy_i for anisotropy_i in anisotropy_coef]

    gauss1 = stack.gaussian_filter(image, corrected_sigma_1)
    gauss2 = stack.gaussian_filter(image, corrected_sigma_2)
    res = _unsigned_int_substract(gauss1, gauss2)

    return res