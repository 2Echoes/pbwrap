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

def correlation_score(pearson_correlation, random_correlation) : 
    return (pearson_correlation - random_correlation) / (pearson_correlation + random_correlation)
    