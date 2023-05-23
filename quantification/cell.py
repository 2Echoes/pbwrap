"""
This submodule contains functions to compute features related to cell wide measurement.
"""

import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
import CustomPandasFramework.operations as dataOp
from scipy.ndimage import distance_transform_edt
from bigfish.stack import mean_projection, maximum_projection, check_parameter, check_array
from .utils import unzip
from bigfish.classification import compute_features, get_features_name
from .measures import count_spots_in_mask, compute_mask_area, compute_signalmetrics
from pbwrap.utils import from_label_get_centeroidscoords


def compute_Cell(acquisition_id, cell, pbody_label, dapi, voxel_size = (300,103,103)):
    """Returns DataFrame with expected Cell datashape containing all cell level features. Features are computed using bigFish built in functions.
    
    Parameters
    ----------
        acquisition_id : int
            Unique identifier for current acquisition.
        cell_id : int 
            Unique identifier for current cell.
        cell : dict
            Dictionary computed from bigFish.multistack.extract_cell
    
    Returns
    -------
        new_Cell : pd.Dataframe
    """
    #Integrity checks
    check_parameter(acquisition_id = (int), cell = (dict), voxel_size = (tuple, list), dapi = (np.ndarray), pbody_label = (np.ndarray))
    check_array(dapi, ndim=3)
    check_array(pbody_label, ndim= 2) # TODO : update if 3D seg is performed for pbody_label.

    #Extracting bigfish cell information
    voxel_size_yx = voxel_size[1]
    cell_mask = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"] 
    rna_coord = cell["rna_coord"]
    foci_coord = cell["foci"]
    ts_coord = cell["transcription_site"]
    malat1_coord = cell["malat1_coord"]
    smfish = cell["smfish"]
    min_y, min_x, max_y, max_x = cell["bbox"]

    #Computing pbody coords from masks
    assert cell_mask.dtype == bool, "cell_mask is not boolean this should NOT happen."
    pbody_label: np.ndarray = pbody_label[min_y : max_y, min_x : max_x]
    pbody_label[~cell_mask] = 0 # Excluding p-bodies in the neighborhood but not in the cell
    pbody_mask = pbody_label.astype(bool)
    centroids_dict = from_label_get_centeroidscoords(pbody_label)
    pbody_centroids = list(zip(centroids_dict["centroid-0"], centroids_dict["centroid-1"]))

    #Removing p-bodies with centroids outside of the cell even though part of the mask is inside cell
    for centroid in pbody_centroids :
        y,x = int(centroid[0]), int(centroid[1]) #pbody centroid is determined with +- 0.5 px error
        if not cell_mask[y,x] :
            pbody_centroids.remove(centroid)
    pbody_centroids = np.array(pbody_centroids,  dtype= int)

    has_pbody = pbody_centroids.ndim > 1
    del centroids_dict 

    #BigFish built in features
    if not has_pbody:
        features, features_names = compute_features(cell_mask= cell_mask, nuc_mask= nuc_mask, ndim= 3, rna_coord= rna_coord, smfish= smfish, foci_coord= foci_coord, voxel_size_yx= voxel_size_yx,
        compute_centrosome=False,
        compute_distance=True,
        compute_intranuclear=True,
        compute_protrusion=True,
        compute_dispersion=True,
        compute_topography=True,
        compute_foci=True,
        compute_area=True,
        return_names=True)
        
    #if there is pbody
    else:
        features, features_names = compute_features(cell_mask= cell_mask, nuc_mask= nuc_mask, ndim= 3, rna_coord= rna_coord, smfish= smfish, centrosome_coord= pbody_centroids, foci_coord= foci_coord, voxel_size_yx= voxel_size_yx,
            compute_centrosome=True,
            compute_distance=True,
            compute_intranuclear=True,
            compute_protrusion=True,
            compute_dispersion=True,
            compute_topography=True,
            compute_foci=True,
            compute_area=True,
            return_names=True)
    
    #Custom features
    cluster_number = len(ts_coord) + len(foci_coord)
    nucleus_area_px = compute_mask_area(nuc_mask, unit= 'px', voxel_size= voxel_size)
    nucleus_area_nm = compute_mask_area(nuc_mask, unit= 'nm', voxel_size= voxel_size)
    #signal features
    nucleus_mip_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mip')
    nucleus_mean_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mean')

    #Adding custom signal features to DataFrame
    features = np.append(features, [nucleus_mip_signal_metrics["mean"], nucleus_mip_signal_metrics["max"], nucleus_mip_signal_metrics["min"], nucleus_mip_signal_metrics["median"],
                                    nucleus_mean_signal_metrics["mean"], nucleus_mean_signal_metrics["max"], nucleus_mean_signal_metrics["min"], nucleus_mean_signal_metrics["median"]])
    
    features_names += ["nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                        "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]

    #malat features
    malat1_spot_in_nuc = count_spots_in_mask(malat1_coord, nuc_mask)
    malat1_spot_in_cyto = count_spots_in_mask(malat1_coord, cell_mask) - malat1_spot_in_nuc
    
    #pbody features
    pbody_num = len(pbody_centroids)
    if has_pbody :
        pbody_area_px = compute_mask_area(pbody_mask, unit= 'px', voxel_size= voxel_size)
        pbody_area_nm = compute_mask_area(pbody_mask, unit= 'nm', voxel_size= voxel_size)
        rna_spot_in_pbody = count_spots_in_mask(rna_coord, pbody_mask)
        count_pbody_nucleus = count_spots_in_mask(pbody_centroids, nuc_mask)
        count_pbody_cytoplasm = count_spots_in_mask(pbody_centroids, cell_mask) - count_pbody_nucleus
        pbody_closer_than_1000_nm = count_rna_close_pbody(pbody_mask= pbody_mask, spots_coords= rna_coord, distance_nm= 1000, voxel_size= voxel_size)
        pbody_closer_than_1500_nm = count_rna_close_pbody(pbody_mask= pbody_mask, spots_coords= rna_coord, distance_nm= 1500, voxel_size= voxel_size)
        pbody_closer_than_2000_nm = count_rna_close_pbody(pbody_mask= pbody_mask, spots_coords= rna_coord, distance_nm= 2000, voxel_size= voxel_size)
    else :
        pbody_area_px = np.NaN
        pbody_area_nm = np.NaN
        rna_spot_in_pbody = np.NaN
        count_pbody_nucleus = np.NaN
        count_pbody_cytoplasm = np.NaN
        pbody_closer_than_1000_nm = np.NaN
        pbody_closer_than_1500_nm = np.NaN
        pbody_closer_than_2000_nm = np.NaN


    #Adding custom features to DataFrames
    features = np.append(features, [malat1_spot_in_nuc, malat1_spot_in_cyto, cluster_number,nucleus_area_px,nucleus_area_nm,
                         rna_spot_in_pbody, pbody_num, pbody_area_px, pbody_area_nm, count_pbody_nucleus, count_pbody_cytoplasm, pbody_closer_than_1000_nm, pbody_closer_than_1500_nm, pbody_closer_than_2000_nm])
    features_names += ['malat1 spots in nucleus', 'malat1 spots in cytoplasm', 'cluster number','nucleus area (px)','nucleus area (nm^2)',
               'rna spots in pbody', 'pbody number', 'pbody area (px)', 'pbody area (nm^2)', "count pbody in nucleus", "count pbody in cytoplasm", "rna 1000 nm from pbody", "rna 1500 nm from pbody", "rna 2000 nm from pbody"]
    header = ["id", "AcquisitionId"] + features_names
    data = np.append([0, acquisition_id], features)
    data = data.reshape([1,-1])

    #Ensuring correct datashape
    datashape_ref = DataFrame.newframe_Cell()
    new_Cell = pd.DataFrame(data= data, columns= header)
    new_Cell["plot index"] = np.NaN
    if not has_pbody :
        for feature in get_features_name(names_features_centrosome= True) :
            new_Cell[feature] = np.NaN
    dataOp.check_samedatashape(new_Cell, datashape_ref) # Ensure datashape stability along different runs

    return new_Cell




def nucleus_signal_metrics(cell, channel, projtype = 'mip') :
    """Returns dict containing signal related measures : 'min', 'max', '1 percentile', '9 percentile', 'mean' and 'median'.
      Computed from channel signal in cell's nucleus mask. Signal measures are computed from 2D cell, so channel is projected z-wise according to projtype (provided channel is 3D).
    
        Parameters
        ----------
            cell : dict
                Dictionary computed from bigFish.multistack.extract_cell
            channel : np.ndarray
                Channel from which intensity is computed
            projtype : str
                can either be 'mip' or 'mean'.

        Returns
        -------
            mean_sig : float
        
    """
    min_y, min_x, max_y, max_x = cell["bbox"]
    channel_cropped = channel[:, min_y:max_y, min_x:max_x]
 

    if channel.ndim == 3 :
        if projtype == 'mip' : 
            channel_crop_proj = maximum_projection(channel_cropped)
        elif projtype == 'mean' :
            channel_crop_proj = mean_projection(channel_cropped)
    
    nucleus_mask = cell["nuc_mask"]
    metrics = compute_signalmetrics(channel_crop_proj, nucleus_mask)
    return metrics



def count_rna_close_pbody(pbody_mask: np.ndarray, spots_coords: 'list[tuple]', distance_nm: float, voxel_size: 'tuple[float]')-> int :
    """
    Count number of RNA (spots) closer than 'distance_nm' from a p-body (mask).
    """
    
    check_parameter(pbody_mask = (np.ndarray), spots_coords = (list, np.ndarray), distance_nm = (int, float), voxel_size = (tuple, list))

    if pbody_mask.ndim != 2: raise ValueError("Unsupported p_body mask dimension. Only 2D arrays are supported.")
    if type(spots_coords) == np.ndarray : spots_coords = list(spots_coords)
    if len(voxel_size) == 3 :
        y_scale = voxel_size[1]
        x_scale = voxel_size[2]
    elif len(voxel_size) == 2 :
        y_scale = voxel_size[0]
        x_scale = voxel_size[1]
    else : raise ValueError("Incorrect voxel_size length should be either 2 or 3. {0} was given".format(len(voxel_size)))


    frompbody_distance_map = distance_transform_edt(np.logical_not(pbody_mask), sampling= [y_scale, x_scale])
    rna_distance_map = np.ones_like(pbody_mask) * -999
    if len(spots_coords) == 0 : return 0
    if len(spots_coords[0]) == 2 :
        y_coords, x_coords = unzip(spots_coords)
    elif len(spots_coords[0]) == 3 :
        z_coords, y_coords, x_coords = unzip(spots_coords)
        del z_coords
    else : 
        z_coords, y_coords, x_coords,*_ = unzip(spots_coords)
        del z_coords,_
    rna_distance_map[y_coords, x_coords] = frompbody_distance_map[y_coords, x_coords] # This distance maps gives the distance of each RNA to the closest p-body
    count_map = rna_distance_map[rna_distance_map >= 0] <= distance_nm
    
    values,count = np.unique(count_map, return_counts= True)
    if not True in values : 
        count = 0
    else:
        index = list(values).index(True)
        count = count[index]
    
    return count

        