import numpy as np
import pandas as pd
import bigfish.classification as classification
import bigfish.plot as plot
from skimage.measure import regionprops_table
from .measures import compute_mask_area, nucleus_signal_metrics
from ._errors import QuantificationError
from ..detection.centrosome import detect_centrosome

def compute_Cell(cell, voxel_size, dapi_stack, acquisition_id, pipeline_name: str ,**pipeline_parameters) :

    if 'centrosome' in pipeline_name :

        if 'centrosome_presegmentation' not in pipeline_parameters : raise KeyError("For centrosome pipeline, 'centrosome_presegmentation' parameter is expected.")
        centrosome_presegmentation = pipeline_parameters['centrosome_presegmentation']
        new_cell = _centrosome_cell_quant(cell, voxel_size, dapi_stack, acquisition_id, centrosome_presegmentation=centrosome_presegmentation)
    
    else : raise ValueError("pipeline name not recognised.")

    return new_cell


def _main_cell_quant(cell, voxel_size, dapi_stack, acquisition_id, compute_centrosome= False, centrosome_coords= None) :
    """
    Basic function for cell quantification using bigfish.

    Parameters
    ----------
        cell : dict  
            computed from `bigfish.multistack.extract_cell`
        voxel_size : tuple(z,y,x)
        dapi_stack : np.ndarray
            3D uncropped dapi stack.

    Returns
    -------
        Cell : pd.DataFrame  
    """
    
    #Extracting bigfish cell information
    voxel_size_yx = voxel_size[1]
    cell_mask: np.ndarray = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"] 
    rna_coord = cell["rna_coord"]
    smfish = cell["smfish"]
    min_y, min_x, max_y, max_x = cell["bbox"]
    label = cell["cell_id"] # is the label of this cell in cell_label

    if "transcription_site" in cell.keys() :
        ts_coord = cell["transcription_site"]
    else : ts_coord = np.empty(shape=(0,0), dtype= np.int64)

    if "foci" in cell.keys() :
        foci_coord = cell["foci"]
    else : foci_coord = np.empty(shape=(0,0), dtype= np.int64)

    if not isinstance(centrosome_coords, (np.ndarray, list)) and compute_centrosome:
        print("Warning : compute_centrosome is set to True but centrosome_coords parameter is not np.ndarray nor list. \nCompute centrosome has been set to False; centrosome_coords type : {0}".format(type(centrosome_coords)))
        compute_centrosome = False
    if type(centrosome_coords) == list : centrosome_coords = np.array(centrosome_coords, dtype=int)
    features, features_names = classification.compute_features(cell_mask= cell_mask, nuc_mask= nuc_mask, ndim= 3, rna_coord= rna_coord, smfish= smfish, foci_coord= foci_coord, voxel_size_yx= voxel_size_yx, centrosome_coord= centrosome_coords,
        compute_centrosome=compute_centrosome,
        compute_distance=True,
        compute_intranuclear=True,
        compute_protrusion=True,
        compute_dispersion=True,
        compute_topography=True,
        compute_foci=True,
        compute_area=True,
        return_names=True
    )

    #Custom features
    cell_props_table = regionprops_table(cell_mask.astype(int), properties= ["centroid"])
    cell_coordinates = (float(cell_props_table["centroid-0"] + min_y), float(cell_props_table["centroid-1"] + min_x))
    del cell_props_table
    cluster_number = len(ts_coord) + len(foci_coord)
    nucleus_area_px = compute_mask_area(nuc_mask, unit= 'px', voxel_size= voxel_size)
    nucleus_area_nm = compute_mask_area(nuc_mask, unit= 'nm', voxel_size= voxel_size)
    #signal features
    nucleus_mip_signal_metrics = nucleus_signal_metrics(cell, channel= dapi_stack, projtype= 'mip')
    nucleus_mean_signal_metrics = nucleus_signal_metrics(cell, channel= dapi_stack, projtype= 'mean')
 
    #Adding custom features to DataFrame
    features = list(features)
    features.extend([cell_coordinates, label, cell["bbox"], cluster_number, nucleus_area_px, nucleus_area_nm,
                         nucleus_mip_signal_metrics["mean"], nucleus_mip_signal_metrics["max"], nucleus_mip_signal_metrics["min"], nucleus_mip_signal_metrics["median"],
                         nucleus_mean_signal_metrics["mean"], nucleus_mean_signal_metrics["max"], nucleus_mean_signal_metrics["min"], nucleus_mean_signal_metrics["median"]])
     
    features_names += [ "cell_coordinates", "label", "bbox", "cluster_number", "nucleus_area_px", "nucleus_area_nm",
                        "nucleus_mip_mean_signal","nucleus_mip_max_signal","nucleus_mip_min_signal","nucleus_mip_median_signal",
                        "nucleus_mean_mean_signal","nucleus_mean_max_signal","nucleus_mean_min_signal","nucleus_mean_median_signal"]
    
    data = [acquisition_id] + features
    header = ['AcquisitionId'] + features_names

    new_cell = pd.DataFrame(data= [data], columns= header)
    return new_cell


def _centrosome_cell_quant(cell, voxel_size, dapi_stack, acquisition_id, centrosome_presegmentation) :
    
    centrosome_coords = detect_centrosome(cell=cell, centrosome_presegmentation= centrosome_presegmentation)
    centrosome_number = len(centrosome_coords)

    if centrosome_number == 0 : raise QuantificationError("No centrosome found")
    else :
        try :
            cell_res = _main_cell_quant(cell=cell, voxel_size=voxel_size, dapi_stack=dapi_stack, acquisition_id=acquisition_id, compute_centrosome= True, centrosome_coords=centrosome_coords)
            cell_res.at[0, "centrosome_number"] = centrosome_number
            centrosome_coords = tuple([tuple(coords) for coords in centrosome_coords])
            cell_res["centrosome_coords"] = np.array(tuple([tuple(coords) for coords in centrosome_coords]),dtype= int)
            print(cell_res)
        except IndexError as error :
            print("centrosome coords : ", centrosome_coords)
            raise error
        return cell_res