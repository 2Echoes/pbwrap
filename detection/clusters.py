"""
Submodule aiming to handle functions around bigfish cluster detection.
"""

import numpy as np
import pandas as pd
import bigfish.detection as detection
from bigfish.stack import check_parameter
from ..utils import nanometer_to_pixel



def _compute_clustered_spots_dataframe(clustered_spots) :
    if len(clustered_spots) == 0 : return pd.DataFrame(columns= ["id", "cluster_id", "z", "y", "x"])
    z, y ,x, cluster_index = list(zip(*clustered_spots))
    ids = np.arange(len(clustered_spots))

    df = pd.DataFrame({
        "id" : ids
        ,"cluster_id" : cluster_index
        ,"z" : z
        ,"y" : y
        ,"x" : x
    })
    null_idx = df[df['cluster_id'] == -1].index
    df.loc[null_idx, 'cluster_id'] = np.NaN

    return df

def _compute_cluster_dataframe(clusters) :
    if len(clusters) == 0 : return pd.DataFrame(columns= ["id", "z", "y", "x", "spot_number"])
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
        res = {'clustered_spots' : [], 'clusters' : [], 'clustered_spots_dataframe' : pd.DataFrame(columns= ["id", "cluster_id", "z", "y", "x"]), 'clusters_dataframe' : ["id", "z", "y", "x", "spot_number"]}
        return {key : res[key] for key in keys_to_compute}
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


def get_centroids_list(clusters_df) :

    """
    clusters_list should be a pd.DataFrame with ['z', 'y', 'x'] or ['y', 'x'] keys.
    """

    if 'y' in clusters_df.columns and 'x' in clusters_df.columns :
        if 'z' in clusters_df.columns : keys = [clusters_df['z'], clusters_df['y'], clusters_df['x']]
        else : keys = [clusters_df['y'], clusters_df['x']]
    else : raise ValueError("Expected keys : ['z', 'y', 'x'] or ['y', 'x']")

    return list(zip(*keys))


def _compute_critical_spot_number(xy_pixel_radius, z_pixel_radius, density) :
    res = 4/3*np.pi*np.square(xy_pixel_radius)*z_pixel_radius*density/100
    return int(round(res))


def remove_artifact(deconvoluted_spots, artifact_radius, voxel_size , spot_density = 20) :
    """
    Artifact are detected as spherical clusters of radius 'artifact_size' and with an average density 'spot_density' of spot within the cluster.
    All spots within the artifact are then removed from deconvoluted_spos.
    
    Critical number of spot is computed as :
    >>> (total_pixel_approximation) * spot_density /100
    >>> with total_pixel_approximation = 4/3*pi*(artifact_radius_xy)²*artifact_radius_z ~rounded to unity

    Parameters
    ----------
        deconvoluted_spots : np.ndarray(z,y,x)
            A dense region decomposition is highly recommended prior to this function
        artifact_radius : int
            in nm
        voxel_size : tuple (z,y,x)
        spot_density : float 
            in range ]0,100]
    """
    
    z_pixel_radius , xy_pixel_radius = nanometer_to_pixel([artifact_radius, artifact_radius], scale= voxel_size[:2])
    critical_spot_number = _compute_critical_spot_number(xy_pixel_radius=xy_pixel_radius, z_pixel_radius=z_pixel_radius, density=spot_density)
    artifacts_df:pd.DataFrame = cluster_detection(deconvoluted_spots, voxel_size=voxel_size, radius= artifact_radius, nb_min_spots=critical_spot_number, keys_to_compute= ['clusters_dataframe'])['clusters_dataframe']
    drop_index = artifacts_df[artifacts_df["cluster_id"].isna()].index
    artifacts_df = artifacts_df.drop(drop_index, axis= 0)

    clean_spots = get_centroids_list(artifacts_df)

    return np.array(clean_spots, dtype= int)