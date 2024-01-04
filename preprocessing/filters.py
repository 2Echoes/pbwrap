import bigfish.stack as stack
import numpy as np

def variance_filter(image, kernel_size= 1, kernel_shape='disk') :
    stack.check_parameter(image, np.ndarray)
    squared_mean_image = np.power(stack.mean_filter(image, kernel_shape=kernel_shape, kernel_size=kernel_size), 2)
    mean_squared_image = stack.mean_filter(np.power(image, 2), kernel_shape=kernel_shape, kernel_size=kernel_size)

    return squared_mean_image - mean_squared_image