import numpy as np
import bigfish.stack as stack
from ..utils import check_parameter

def pearson_correlation(image1: np.ndarray, image2: np.ndarray) :

    """
    Returns pearson correlation coefficient between image 1 and image 2.
    """

    check_parameter(image1= np.ndarray, image2= np.ndarray)

    X1 = image1.flatten()
    X2 = image2.flatten()

    pearson_matrix = np.corrcoef(X1,X2)

    return pearson_matrix[0,1]

def detection_correlation(signal, spots) : 
    if len(spots) == 0 : return 0
    reconstructed_signal = np.zeros_like(signal)
    