import numpy as np
from skimage.segmentation import find_boundaries
from skimage.measure import regionprops_table

def compute_Pbody(AcquisitionId,Pbody_label) :
    """
    
    """
    Pbody_dictionary = compute_Pbody_dictionary(Pbody_label)
    nbre_pbody = len(Pbody_dictionary["id"])
    ids = np.arange(nbre_pbody)
    AcquisitionIds = [AcquisitionId] * nbre_pbody
    CellIds = np.nan # Computed in 'insert update function name'

    pass



def compute_Pbody_dictionary(Pbody_label):

    Pbody_dictionary = regionprops_table(Pbody_label, properties= ['label', 'centroid','area'])
    boundary_coordinates_list = []
    for label in Pbody_dictionary['label'] :
        print(label)
        mask = np.zeros_like(Pbody_label)
        mask[Pbody_label == label] = 1
        boundary_coordinates_list += [boundary_coordinates(mask)]
    Pbody_dictionary["boudary_coordinates"] = boundary_coordinates_list
    return Pbody_dictionary


def boundary_coordinates(mask, to_tuple = True) :

    boundary = find_boundaries(mask, background= 0, mode = 'inner')
    boundary_coordinates = np.argwhere(boundary).tolist()

    if to_tuple: return tuple(map(tuple, boundary_coordinates))
    else : return boundary_coordinates
