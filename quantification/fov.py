"""
This submodule contains functions to compute features related to fov wide measurement.
"""
import numpy as np
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
from bigfish.stack import mean_projection
from .measures import compute_signalmetrics


def compute_fov(acquisition_id, rootfilename:str, rna_name:str, cell_number: int, rna_threshold:float, malat_threshold:float, dapi: np.ndarray, nucleus_mask: np.ndarray) :
    """Returns DataFrame with expected Acquisition datashape containing all acquisition (fov) level features."""

    if dapi.ndim == 3 : 
        dapi = mean_projection(dapi)
    elif dapi.ndim == 2 :
        pass
    else: raise ValueError("dapi signal should either be a 2D or 3D np.ndarray.")
    if nucleus_mask.dtype != bool : raise TypeError("nucleus mask should be of dtype bool.")
    
    
    dapi_backgrnd_metrics = compute_signalmetrics(dapi, ~nucleus_mask)

    new_Acquisition = DataFrame.newframe_Acquisitions()
    new_Acquisition["id"] = [acquisition_id]
    new_Acquisition["rootfilename"] = [rootfilename]
    new_Acquisition["rna name"] = [rna_name]
    new_Acquisition["cell number"] = [cell_number]
    new_Acquisition["RNA spot threshold"] = [rna_threshold]
    new_Acquisition["malat1 spot threshold"] = [malat_threshold]
    new_Acquisition["background dapi signal mean"] = dapi_backgrnd_metrics['mean']
    new_Acquisition["background dapi signal std"] = dapi_backgrnd_metrics['std']
    
    
    return new_Acquisition