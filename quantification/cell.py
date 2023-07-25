"""
This submodule contains functions to compute features related to cell wide measurement.
"""
import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
from CustomPandasFramework.integrity import check_samedatashape
from scipy.ndimage import distance_transform_edt
from skimage.measure import regionprops_table
from bigfish.stack import mean_projection, maximum_projection, check_parameter, check_array
from .utils import unzip
from bigfish.classification import compute_features, get_features_name
from .measures import count_spots_in_mask, compute_mask_area, compute_signalmetrics
from pbwrap.utils import from_label_get_centeroidscoords
"""

def compute_Cell(acquisition_id, cell, pbody_label, dapi, voxel_size = (300,103,103)):
    """
"""
    Returns DataFrame with expected Cell datashape containing all cell level features. Features are computed using bigFish built in functions.
    
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
"""
    #Integrity checks
    check_parameter(acquisition_id = (int), cell = (dict), voxel_size = (tuple, list), dapi = (np.ndarray), pbody_label = (np.ndarray))
    check_array(dapi, ndim=3)
    check_array(pbody_label, ndim= 2) # TODO : update if 3D seg is performed for pbody_label.

    #Extracting bigfish cell information
    voxel_size_yx = voxel_size[1]
    cell_mask: np.ndarray = cell["cell_mask"]
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
    assert pbody_label.shape == cell_mask.shape
    pbody_copy = pbody_label.copy()
    pbody_copy[~cell_mask] = 0 # Excluding p-bodies in the neighborhood but not in the cell
    pbody_mask = pbody_copy.astype(bool)
    centroids_dict = from_label_get_centeroidscoords(pbody_label)

    Y,X = centroids_dict["centroid-0"], centroids_dict["centroid-1"]
    Y,X = np.array(Y).round().astype(int) , np.array(X).round().astype(int)

    Y_array = Y[cell_mask[Y,X]]
    X_array = X[cell_mask[Y,X]]

    Y_abs,X_abs = Y_array + min_y, X_array + min_x 
    pbody_coordinates = list(zip(Y_abs, X_abs))
    
    pbody_centroids = np.array(list(zip(Y_array, X_array)))
    pbody_num = count_spots_in_mask(pbody_centroids, cell_mask)
    has_pbody = pbody_num > 0
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
    features = list(features)
    
    #Custom features
    cell_props_table = regionprops_table(cell_mask.astype(int), properties= ["centroid"])
    cell_coordinates = (float(cell_props_table["centroid-0"] + min_y), float(cell_props_table["centroid-1"] + min_x))
    label = cell["cell_id"] # is the label of this cell in cell_label
    del cell_props_table
    cluster_number = len(ts_coord) + len(foci_coord)
    nucleus_area_px = compute_mask_area(nuc_mask, unit= 'px', voxel_size= voxel_size)
    nucleus_area_nm = compute_mask_area(nuc_mask, unit= 'nm', voxel_size= voxel_size)
    #signal features
    nucleus_mip_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mip')
    nucleus_mean_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mean')

    #Adding custom signal features to DataFrame
    features.extend([cell_coordinates, label, cell["bbox"], pbody_coordinates,
                         nucleus_mip_signal_metrics["mean"], nucleus_mip_signal_metrics["max"], nucleus_mip_signal_metrics["min"], nucleus_mip_signal_metrics["median"],
                         nucleus_mean_signal_metrics["mean"], nucleus_mean_signal_metrics["max"], nucleus_mean_signal_metrics["min"], nucleus_mean_signal_metrics["median"]])
    
    features_names += [ "cell_coordinates", "label", "bbox", "pbody coordinates",
                        "nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                        "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]

    #malat features
    malat1_spot_in_nuc = count_spots_in_mask(malat1_coord, nuc_mask)
    malat1_spot_in_cyto = count_spots_in_mask(malat1_coord, cell_mask) - malat1_spot_in_nuc
    
    #pbody features
    if has_pbody :
        pbody_area_px = compute_mask_area(pbody_mask, unit= 'px', voxel_size= voxel_size)
        pbody_area_nm = compute_mask_area(pbody_mask, unit= 'nm', voxel_size= voxel_size)
        rna_spot_in_pbody = count_spots_in_mask(rna_coord, pbody_mask)
        count_pbody_nucleus = count_spots_in_mask(pbody_centroids, nuc_mask)
        count_pbody_cytoplasm = pbody_num - count_pbody_nucleus
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
    features.extend([malat1_spot_in_nuc, malat1_spot_in_cyto, cluster_number,nucleus_area_px,nucleus_area_nm,
                         rna_spot_in_pbody, pbody_num, pbody_area_px, pbody_area_nm, count_pbody_nucleus, count_pbody_cytoplasm, pbody_closer_than_1000_nm, pbody_closer_than_1500_nm, pbody_closer_than_2000_nm])
    features_names += ['malat1 spots in nucleus', 'malat1 spots in cytoplasm', 'cluster number','nucleus area (px)','nucleus area (nm^2)',
               'rna spots in pbody', 'pbody number', 'pbody area (px)', 'pbody area (nm^2)', "count pbody in nucleus", "count pbody in cytoplasm", "rna 1000 nm from pbody", "rna 1500 nm from pbody", "rna 2000 nm from pbody"]
    header = ["id", "AcquisitionId"] + features_names
    data = [0, acquisition_id] 
    data.extend(features)

    #Ensuring correct datashape
    datashape_ref = DataFrame.newframe_Cell()
    new_Cell = pd.DataFrame(data= [data], columns= header)
    new_Cell["plot index"] = np.NaN
    if not has_pbody :
        for feature in get_features_name(names_features_centrosome= True) :
            new_Cell[feature] = np.NaN
    check_samedatashape(new_Cell, datashape_ref) # Ensure datashape stability along different runs
    return new_Cell

"""
def compute_Nucleus(cell: dict, dapi, voxel_size, acquisition_id) : 
    """
    Is a subpart of compute Cell. Meaning that it contains the fraction of measure of 'compute_Cell' which are related to DAPI signal (Nucleus).
    When calling multistack.extract_cell use nucleus mask as cell mask.
    """
    

    new_Cell_columns = DataFrame.newframe_Cell(names_features_distance= False,
                    names_features_area= False,
                    names_features_centrosome= False,
                    names_features_dispersion= False,
                    names_features_foci= False,
                    names_features_intranuclear= False,
                    names_features_protrusion= False,
                    names_features_topography= False,
                    names_features_signal= True,
                    names_features_malat1= False,
                    names_features_pbody= False,
                    plot_index= False).columns
    
    min_y, min_x, max_y, max_x = cell["bbox"]
    nuc_mask = cell["cell_mask"]
    cell_props_table = regionprops_table(nuc_mask.astype(int), properties= ["centroid"])
    cell_coordinates = (float(cell_props_table["centroid-0"] + min_y), float(cell_props_table["centroid-1"] + min_x))

    label = cell["cell_id"] # is the label of this cell in cell_label
    del cell_props_table
    nucleus_area_px = compute_mask_area(nuc_mask, unit= 'px', voxel_size= voxel_size)
    nucleus_area_nm = compute_mask_area(nuc_mask, unit= 'nm', voxel_size= voxel_size)
    #signal features
    nucleus_mip_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mip', use_cell_mask= True)
    nucleus_mean_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mean', use_cell_mask= True)

    features = [cell_coordinates, label, cell["bbox"],
                         nucleus_mip_signal_metrics["mean"], nucleus_mip_signal_metrics["max"], nucleus_mip_signal_metrics["min"], nucleus_mip_signal_metrics["median"],
                         nucleus_mean_signal_metrics["mean"], nucleus_mean_signal_metrics["max"], nucleus_mean_signal_metrics["min"], nucleus_mean_signal_metrics["median"],
                         nucleus_area_px, nucleus_area_nm]


    data = [0, acquisition_id] 
    data.extend(features)
    new_Cell = pd.DataFrame(columns= new_Cell_columns, data= [data])

    return new_Cell

def nucleus_signal_metrics(cell, channel, projtype = 'mip', use_cell_mask= False) :
    """
    Returns dict containing signal related measures : 'min', 'max', '1 percentile', '9 percentile', 'mean' and 'median'.
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
    
    if use_cell_mask : nucleus_mask = cell["cell_mask"]
    else : nucleus_mask = cell["nuc_mask"]

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





########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################



"""
This submodule contains functions to compute features related to cell wide measurement.
"""
import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
from CustomPandasFramework.integrity import check_samedatashape
from scipy.ndimage import distance_transform_edt
from skimage.measure import regionprops_table
from bigfish.stack import mean_projection, maximum_projection, check_parameter, check_array
from .utils import unzip
from bigfish.classification import compute_features, get_features_name
from .measures import count_spots_in_mask, compute_mask_area, compute_signalmetrics
from pbwrap.utils import from_label_get_centeroidscoords
# 
# 
def compute_Cell(acquisition_id, cell, Pbody_Acquisition:pd.DataFrame, dapi, cell_label, voxel_size = (300,103,103)):
    """
    Returns DataFrame with expected Cell datashape containing all cell level features. Features are computed using bigFish built in functions.
    # 
    Parameters
    ----------
        acquisition_id : int
            Unique identifier for current acquisition.
        cell_id : int 
            Unique identifier for current cell.
        cell : dict
            Dictionary computed from bigFish.multistack.extract_cell
    # 
    Returns
    -------
        new_Cell : pd.Dataframe
    """
    #Integrity checks
    check_parameter(acquisition_id = (int), cell = (dict), voxel_size = (tuple, list), dapi = (np.ndarray), Pbody_Acquisition = (pd.DataFrame))
    check_array(dapi, ndim=3)
# 


    #Extracting bigfish cell information
    voxel_size_yx = voxel_size[1]
    cell_mask: np.ndarray = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"] 
    rna_coord = cell["rna_coord"]
    foci_coord = cell["foci"]
    ts_coord = cell["transcription_site"]
    malat1_coord = cell["malat1_coord"]
    smfish = cell["smfish"]
    min_y, min_x, max_y, max_x = cell["bbox"]
    label = cell["cell_id"] # is the label of this cell in cell_label
    
    if Pbody_Acquisition.query('cell_label == {0}'.format(label)).loc[:,"centroid_coordinates"].empty :
        Y_abs, X_abs = np.array([]), np.array([])
        Y, X = np.array([]), np.array([])

    else : 
        Y_abs, X_abs = zip(*Pbody_Acquisition.query('cell_label == {0}'.format(label)).loc[:,"centroid_coordinates"])
        assert len(Y_abs) == len(X_abs)
        Y,X = np.array(Y_abs) - min_y, np.array(X_abs) - min_x
    #Computing pbody coords from masks
    
    assert cell_mask.dtype == bool, "cell_mask is not boolean this should NOT happen."

    pbody_coordinates = list(zip(Y_abs, X_abs))    
    pbody_centroids = np.array(list(zip(Y, X)))
    pbody_num = count_spots_in_mask(pbody_centroids, cell_mask)
    has_pbody = pbody_num > 0

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
    features = list(features)
     
    #Custom features
    cell_props_table = regionprops_table(cell_mask.astype(int), properties= ["centroid"])
    cell_coordinates = (float(cell_props_table["centroid-0"] + min_y), float(cell_props_table["centroid-1"] + min_x))
    label_bis = cell_label[int(cell_coordinates[0]), int(cell_coordinates[1])]
    del cell_props_table
    cluster_number = len(ts_coord) + len(foci_coord)
    nucleus_area_px = compute_mask_area(nuc_mask, unit= 'px', voxel_size= voxel_size)
    nucleus_area_nm = compute_mask_area(nuc_mask, unit= 'nm', voxel_size= voxel_size)
    #signal features
    nucleus_mip_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mip')
    nucleus_mean_signal_metrics = nucleus_signal_metrics(cell, channel= dapi, projtype= 'mean')
 
    #Adding custom signal features to DataFrame
    features.extend([cell_coordinates, label, label_bis, cell["bbox"], pbody_coordinates,
                         nucleus_mip_signal_metrics["mean"], nucleus_mip_signal_metrics["max"], nucleus_mip_signal_metrics["min"], nucleus_mip_signal_metrics["median"],
                         nucleus_mean_signal_metrics["mean"], nucleus_mean_signal_metrics["max"], nucleus_mean_signal_metrics["min"], nucleus_mean_signal_metrics["median"]])
     
    features_names += [ "cell_coordinates", "label","label_bis", "bbox", "pbody coordinates",
                        "nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                        "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]
 
    #malat features
    malat1_spot_in_nuc = count_spots_in_mask(malat1_coord, nuc_mask)
    malat1_spot_in_cyto = count_spots_in_mask(malat1_coord, cell_mask) - malat1_spot_in_nuc
     
    #pbody features
    if has_pbody :
        count_pbody_nucleus = count_spots_in_mask(pbody_centroids, nuc_mask)
        count_pbody_cytoplasm = pbody_num - count_pbody_nucleus

    else :
        count_pbody_nucleus = np.NaN
        count_pbody_cytoplasm = np.NaN



    #Adding custom features to DataFrames
    features.extend([malat1_spot_in_nuc, malat1_spot_in_cyto, cluster_number,nucleus_area_px,nucleus_area_nm,
                         pbody_num, count_pbody_nucleus, count_pbody_cytoplasm])
    features_names += ['malat1 spots in nucleus', 'malat1 spots in cytoplasm', 'cluster number','nucleus area (px)','nucleus area (nm^2)',
               'pbody number', "count pbody in nucleus", "count pbody in cytoplasm"]
    header = ["id", "AcquisitionId"] + features_names
    data = [0, acquisition_id] 
    data.extend(features)

    #Ensuring correct datashape
    datashape_ref = DataFrame.newframe_Cell()
    new_Cell = pd.DataFrame(data= [data], columns= header)
    new_Cell["plot index"] = np.NaN
    if not has_pbody :
        for feature in get_features_name(names_features_centrosome= True) :
            new_Cell[feature] = np.NaN
    check_samedatashape(new_Cell, datashape_ref) # Ensure datashape stability along different runs
    return new_Cell