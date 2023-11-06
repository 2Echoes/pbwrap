"""
Submodule aiming to handle functions around bigfish cluster detection.
"""

import numpy as np
import pandas as pd
import bigfish.detection as detection
from bigfish.stack import check_parameter




def _compute_clustered_spots_dataframe(clustered_spots) :

    z, y ,x, cluster_index = list(zip(*clustered_spots))
    ids = np.arange(len(clustered_spots))

    df = pd.DataFrame({
        "id" : ids
        ,"cluster_id" : cluster_index
        ,"z" : z
        ,"y" : y
        ,"x" : x
    })

    return df

def _compute_cluster_dataframe(clusters) :
    z, y, x, spots_number, cluster_index = list(zip(*clusters))

    df = pd.DataFrame({
        "id" : cluster_index
        ,"z" : z
        ,"y" : y
        ,"x" : x
        ,"spot_number" : spots_number
    })

    return df


def cluster_detection(spots, voxel_size, radius = 350, nb_min_spots = 4, keys_to_compute = ["clustered_spots", "clusters"]) :
    """
    Performs `bigfish.detection.cluster_detection()` to detect clusters.
    Then offers possibility to get results sorted in pandas dataframe.

    Parameters
    ----------
        spots : np.ndarray
            Coordinates of the detected spots with shape (nb_spots, 3) or (nb_spots, 2).
        voxel_size : int, float, Tuple(int, float) or List(int, float)
            Size of a voxel, in nanometer. One value per spatial dimension (zyx or yx dimensions). If it's a scalar, the same value is applied to every dimensions.
        radius : int
            The maximum distance between two samples for one to be considered as in the neighborhood of the other. Radius expressed in nanometer.
        nb_min_spots : int
            The number of spots in a neighborhood for a point to be considered as a core point (from which a cluster is expanded). This includes the point itself.
        keys_to_compute : list[str], str
            keys from (clustered_spots, clusters, clustered_spots_dataframe, clusters_dataframe)
                --> clustered_spots : np.ndarray
                    Coordinates of the detected spots with shape (nb_spots, 4) or (nb_spots, 3). One coordinate per dimension (zyx or yx coordinates) plus the index of the cluster assigned to the spot. If no cluster was assigned, value is -1.
                --> clusters : np.ndarray
                    Array with shape (nb_clusters, 5) or (nb_clusters, 4). One coordinate per dimension for the clusters centroid (zyx or yx coordinates), the number of spots detected in the clusters and its index.
                --> clustered_spots_dataframe
                --> clusters_dataframe
    
    Returns
    -------
        res : dict
            keys : keys from `keys_to_compute` argument : (clustered_spots, clusters, clustered_spots_dataframe, clusters_dataframe)    
    """

    if isinstance(keys_to_compute, str) : keys_to_compute = [keys_to_compute]
    elif isinstance(keys_to_compute, list) : pass
    else : raise TypeError("Wrong type for keys_to_compute. Should be list[str] or str. It is {0}".format(type(keys_to_compute)))
    if len(spots) == 0 :
        return {key : [] for key in keys_to_compute}
    else : res = {}
    clustered_spots, clusters = detection.detect_clusters(spots, voxel_size= voxel_size, radius= radius, nb_min_spots= nb_min_spots)


    if 'clustered_spots' in keys_to_compute :
        res['clustered_spots'] = clustered_spots
        
    if 'clusters' in keys_to_compute : 
        res['clusters'] = clusters

    if 'clustered_spots_dataframe' in keys_to_compute :
        res['clustered_spots_dataframe'] = _compute_clustered_spots_dataframe(clustered_spots)
    
    if 'clusters_dataframe' in keys_to_compute :
        res['clusters_dataframe'] = _compute_cluster_dataframe(clusters)

    return res

def get_centroids_list(clusters_array: np.ndarray) :
    """
    clusters_list should be a np.ndarray with shape = (clusters_number, 5 or 4) depending if 3D or 2D clusters : z,y,x, spot_number_in_cluster,cluster_index.
    """

    if not isinstance(clusters_array, np.ndarray) : raise TypeError("clusters_array should be of type np.ndarray, it is {0}".format(type(clusters_array)))
    if clusters_array.shape[1] == 5 : dim = 3
    elif clusters_array.shape[1] == 4 : dim = 2
    else : raise ValueError("Expected shape for clusters_array is clusters_number, 5 or 4) depending if 3D or 2D clusters : z,y,x, spot_number_in_cluster,cluster_index. It is {0}".format(clusters_array.shape))

    if dim == 3 :
        z, y, x, spot_number, index = zip(*clusters_array)
        return list(zip(z,y,x))
    else :
        y, x, spot_number, index = zip(*clusters_array)
        return list(zip(y,x))