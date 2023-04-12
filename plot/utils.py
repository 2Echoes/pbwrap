import numpy as np
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter


def from_label_get_centeroidscoords(label):
    """Returns dict{"label", "centroid"}"""

    check_parameter(label = (np.ndarray))
    dim = label.ndim
    properties_dic = regionprops_table(label, properties= ["label","centroid"])
    centroid = properties_dic
    return centroid