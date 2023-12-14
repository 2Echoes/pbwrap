"""
This submodule contains functions to compute features related to cell wide measurement.
"""
import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
import bigfish.multistack as multistack
from CustomPandasFramework.integrity import check_samedatashape
from scipy.ndimage import distance_transform_edt
from skimage.measure import regionprops_table
from bigfish.stack import mean_projection, maximum_projection, check_parameter, check_array
from .utils import unzip
from bigfish.classification import compute_features, get_features_name
from .measures import count_spots_in_mask, compute_mask_area, nucleus_signal_metrics, compute_signalmetrics

"""

def compute_Cell(acquisition_id, cell, pbody_label, dapi, voxel_size = (300,103,103)):
    
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
#######################################################################################################################
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
def compute_Cell(acquisition_id, cell, dapi, cell_label, Pbody_Acquisition:pd.DataFrame=None, voxel_size = (300,103,103)):
    """
    Returns DataFrame with expected Cell datashape containing all cell level features. 
    Features are computed using bigFish built in functions.
    

    Parameters
    ----------
    acquisition_id : int
        Unique identifier for current acquisition.
    cell : dict
        Dictionary computed from bigFish.multistack.extract_cell
    dapi :  np.ndarray(ndim=3)
        raw dapi signal.
    cell_label : np.ndarray(ndim= 2, dtype= int)
        complete fov cell labeling.
    Pbody_Acquisition : pd.DataFrame
        for Pbody pipeline. Computed with quant.compute_Pbody
    voxel_size : tuple (z,y,x)
        microscope scale factors
    
    
    Returns
    -------
        new_Cell : pd.Dataframe
    """
    #Integrity checks
    check_parameter(acquisition_id = (int), cell = (dict), voxel_size = (tuple, list), dapi = (np.ndarray), Pbody_Acquisition = (pd.DataFrame, type(None)))
    check_array(dapi, ndim=3)
# 


    #Extracting bigfish cell information
    voxel_size_yx = voxel_size[1]
    cell_mask: np.ndarray = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"] 
    rna_coord = cell["rna_coord"]

    if "transcription_site" in cell.keys() :
        ts_coord = cell["transcription_site"]
    else : ts_coord = np.empty(shape=(0,0), dtype= np.int64)

    if "foci" in cell.keys() :
        foci_coord = cell["foci"]
    else : foci_coord = np.empty(shape=(0,0), dtype= np.int64)

    if "malat1_coord" in cell.keys() :
        malat1_coord = cell["malat1_coord"]
    else : malat1_coord = np.empty(shape=(0,0), dtype= np.int64)

    if "centrosome_coord" in cell.keys() :
        centrosome_coord = cell["centrosome_coord"]
        has_centrosome = True
    else : centrosome_coord = None; has_centrosome = False

    smfish = cell["smfish"]
    min_y, min_x, max_y, max_x = cell["bbox"]
    label = cell["cell_id"] # is the label of this cell in cell_label

    if Pbody_Acquisition.query('cell_label == {0}'.format(label)).loc[:,"centroid_coordinates"].empty or type(Pbody_Acquisition) == type(None):
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
    assert not has_pbody or not has_centrosome, "Compute cell doesn't support both p-bodies and centrosome input."

    #BigFish built in features
    if not has_pbody:
        features, features_names = compute_features(cell_mask= cell_mask, nuc_mask= nuc_mask, ndim= 3, rna_coord= rna_coord, smfish= smfish, foci_coord= foci_coord, centrosome_coord=centrosome_coord, voxel_size_yx= voxel_size_yx,
        compute_centrosome=has_centrosome,
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
    if not has_pbody  and not has_centrosome:
        for feature in get_features_name(names_features_centrosome= True) :
            new_Cell[feature] = np.NaN
    check_samedatashape(new_Cell, datashape_ref) # Ensure datashape stability along different runs
    return new_Cell


def extract_cell(cell_label,
        ndim,
        nuc_label=None,
        rna_coord=None,
        others_coord: dict=None,
        image=None,
        others_image=None,
        remove_cropped_cell=True,
        check_nuc_in_cell=True):
    
    """
    wrapper : remove from other coords empty spot list

    Extract cell-level results for an image.

    The function gathers different segmentation and detection results obtained
    at the image level and assigns each of them to the individual cells.

    Parameters
    ----------
    cell_label : np.ndarray, np.uint or np.int
        Image with labelled cells and shape (y, x).
    ndim : int
        Number of spatial dimensions to consider (2 or 3).
    nuc_label : np.ndarray, np.uint or np.int
        Image with labelled nuclei and shape (y, x). If None, individual
        nuclei are not assigned to each cell.
    rna_coord : np.ndarray
        Coordinates of the detected RNAs with zyx or yx coordinates in the
        first 3 or 2 columns. If None, RNAs are not assigned to individual
        cells.
    others_coord : Dict[np.ndarray]
        Dictionary of coordinates arrays. For each array of the dictionary,
        the different elements are assigned to individual cells. Arrays should
        be organized the same way than spots: zyx or yx coordinates in the
        first 3 or 2 columns, np.int64 dtype, one element per row. Can be used
        to assign different detected elements to the segmented cells along with
        the spots. If None, no others elements are assigned to the individual
        cells.
    image : np.ndarray, np.uint
        Image in 2-d. If None, image of the individual cells are not extracted.
    others_image : Dict[np.ndarray]
        Dictionary of images to crop. If None, no others image of the
        individual cells are extracted.
    remove_cropped_cell : bool
        Remove cells cropped by the FoV frame.
    check_nuc_in_cell : bool
        Check that each nucleus is entirely localized within a cell.

    Returns
    -------
    fov_results : List[Dict]
        List of dictionaries, one per cell segmented in the image. Each
        dictionary includes information about the cell (image, masks,
        coordinates arrays). Minimal information are:

        * `cell_id`: Unique id of the cell.
        * `bbox`: bounding box coordinates with the order (`min_y`, `min_x`,
          `max_y`, `max_x`).
        * `cell_coord`: boundary coordinates of the cell.
        * `cell_mask`: mask of the cell.

    """

    if type(others_coord) != type(None) :
        others_coord_copy = others_coord.copy()
    
        for key, val in others_coord.items() :
            if len(val) == 0 : del others_coord_copy[key]
            elif key == 'malat1_coord' and len(val[0]) == 0 : del others_coord_copy[key]
    else : others_coord_copy = None


    fov_results = multistack.extract_cell(cell_label= cell_label, 
                                                      ndim=ndim, nuc_label= nuc_label, rna_coord= rna_coord, others_coord= others_coord_copy, image= image, others_image= others_image, remove_cropped_cell=remove_cropped_cell, check_nuc_in_cell=check_nuc_in_cell)
    
    return fov_results