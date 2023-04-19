import pandas as pd
import datetime as dt
import re, inspect
import numpy as np
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
import CustomPandasFramework.operations as dataOp
from bigfish.stack import check_parameter,read_image, check_array
from bigfish.classification import compute_features, get_features_name
from ..quantification import mean_signal, count_spots_in_mask
from pbwrap.utils import from_label_get_centeroidscoords




def get_images(path_input, Input, acquisition_index, channels_list= None) :
    """ Open images from an acquisition within the Input DataFrame using bigfish open method. 
        Outputs as a list of image which can be 2D or 3D depending on the file placed in the input folder.
    
    Parameters
    ----------

        path_input : str
            Full path to the input folder.
        Input : pd.DataFrame
            Input DataFrame resulting from get_Input
        acquisition_index : int
            index refering to which acquisition in the Input frame is to get.
        channels_list : List[str]
            List of strings indicating which channels are to be opened. 
            Note that the output list order will correspond to the order of channels_list.
            if None all channels will be output and sorted in alphabetical order A -> Z.

                
    Returns
    -------
    
        images : List[np.ndarray]
        List of images open using bigfish.stack.read_image()
    """

    #Integrity checks
    check_parameter(path_input = (str), Input = (pd.DataFrame), acquisition_index = (int), channels_list = (list, type(None)))

    for string in channels_list : check_parameter(string = (str))

    if channels_list == None :
         channels_list = Input.values_count(subset= "channel").index.tolist()

    images = []
    for channel in channels_list :
         filename = Input.loc[Input["channel"] == channel].reset_index(drop= True).at[acquisition_index, "filename"]
         images += [read_image(path= path_input + filename)]
    return images

def get_acquisition_num(Input) :
    """ Returns the number of acquistion that will be computed from data placed in the input folder.

    Parameters
    ----------

        Input : pd.DataFrame
            Dataframe containing informations on files placed in the input folder. Should be obtained with get_Input()

    Returns
    -------

        acquisition_num : int
            Number of acquisitions that will be computed from files placed in the input folder.
            This number equals last acquisition idx + 1 since in python indexes start at 0.

    """
    #Integrity :
    check_parameter(Input = (pd.DataFrame))

    acquisition_num = Input["acquisition index"].max() + 1
    return acquisition_num


def get_rootfilename(acquisition_index, Input_frame):
    """Returns root filename of an acquisition from Input frame
    
    Parameters
    ----------
        acquisition_index : int
        Input_frame : pd.DataFrame
            Input dataframes are computed from newFrame_Input
    Returns
    -------
        res : str.
    """
    check_parameter(acquisition_index = (int), Input_frame = pd.DataFrame)

    res = Input_frame.value_counts(subset=["root filename", "acquisition index"]).reset_index(drop= False).at[acquisition_index, "root filename"]
    return res

def get_rnaname(acquisition_index, Input_frame):
    """Returns the RNA name of an acquisition from the Input frame

    Parameters
    ----------
        acquisition_index : int
        Input_frame : pd.DataFrame
            Input dataframes are computed from newFrame_Input
    Returns
    -------
        res : str.
    """
    root = get_rootfilename(acquisition_index, Input_frame)
    regex = "(\w*)f\d{2}--"
    res = re.findall(regex,root)[0]
    return res




def get_Cell(acquisition_id, cell, pbody_label, dapi, voxel_size = (300,103,103)):
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
    check_array(pbody_label, ndim= 2) # TODO : code to update if 3D seg is performed for pbody_label.

    #Extracting bigfish cell information
    voxel_size_yx = voxel_size[1] # la résolution dans le plan devrait être toujours x/y indep
    cell_mask = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"]
    rna_coord = cell["rna_coord"]
    foci_coord = cell["foci"]
    malat1_coord = cell["malat1_coord"]
    smfish = cell["smfish"]
    min_y, min_x, max_y, max_x = cell["bbox"]

    #Computing pbody coords from masks
    pbody_label = pbody_label[min_y : max_y, min_x : max_x]
    pbody_mask = np.zeros_like(pbody_label)
    pbody_mask[pbody_label > 0] = 1
    pbody_mask = np.array(pbody_mask, dtype= bool)
    centroids_dict = from_label_get_centeroidscoords(pbody_label)
    pbody_centroids = np.array(list(zip(centroids_dict["centroid-0"], centroids_dict["centroid-1"])), dtype= int)
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
    #signal features
    mip_mean_intensity = mean_signal(cell, channel= dapi, projtype= 'mip')
    mean_mean_intensity = mean_signal(cell, channel= dapi, projtype= 'mean')
    #malat features
    malat1_spot_in_nuc = count_spots_in_mask(malat1_coord, nuc_mask)
    malat1_spot_in_cyto = count_spots_in_mask(malat1_coord, cell_mask) - malat1_spot_in_nuc
    #pbody features
    rna_spot_in_pbody = count_spots_in_mask(rna_coord, pbody_mask)
    pbody_num = len(pbody_centroids)


    #Adding custom features to DataFrames
    features = np.append(features, [mip_mean_intensity, mean_mean_intensity, malat1_spot_in_nuc, malat1_spot_in_cyto, rna_spot_in_pbody, pbody_num])
    features_names += ['Mean Intensity (MIP)', 'Mean Intensity (MeanProj)', 'malat1 spots in nucleus', 'malat1 spots in cytoplasm', 'rna spots in body', 'pbody number']
    header = ["id", "AcquisitionId"] + features_names
    data = np.append([0, acquisition_id], features)
    data = data.reshape([1,-1])

    #Ensuring correct datashape
    datashape_ref = DataFrame.newframe_Cell()
    new_Cell = pd.DataFrame(data= data, columns= header)
    if not has_pbody :
        for feature in get_features_name(names_features_centrosome= True) :
            new_Cell[feature] = np.NaN
    dataOp.check_samedatashape(new_Cell, datashape_ref) # Ensure datashape stability along different runs

    return new_Cell


def from_Acquisition_get_rna(Acquisition: pd.DataFrame) -> pd.Index :
    return list(Acquisition.value_counts(subset= "rna name").index)



def from_rna_get_Cells(rna: 'list[str]', Cell: pd.DataFrame, Acquisition: pd.DataFrame) -> pd.DataFrame :
    """
    Returns sub-table from Cell containing only cells which rna name (from Acquisition) matches one in 'rna'.
    Also 'rna name' column is added to returned Cell sub-table.
    """


    if type(rna) == str : rna = [rna]

    join_frame = dataOp.keep_columns(Dataframe= pd.merge(Cell, Acquisition, how= 'left', left_on= 'AcquisitionId', right_on= 'id'),
                                     columns= ["rna name"] + list(Cell.columns))
    print(join_frame)
    drop_index = join_frame.query('`rna name` not in {0}'.format(rna)).index
    join_frame = join_frame.drop(axis= 0, index= drop_index)

    return join_frame



def _get_varname(var):
    #To be used  within a function.
    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var][0]



def get_datetime():
    return dt.datetime.now().strftime("%Y%m%d %H-%M-%S")