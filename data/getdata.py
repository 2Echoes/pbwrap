import pandas as pd
import datetime as dt
import re, inspect
import numpy as np
import CustomPandasFramework.PBody_project.DataFrames as DataFrame
import CustomPandasFramework.operations as dataOp
from bigfish.stack import check_parameter,read_image
from bigfish.classification import compute_features





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




def get_Cell(acquisition_id, cell, voxel_size = (300,103,103)):
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
    check_parameter(acquisition_id = (int), cell = (dict), voxel_size = (tuple, list))

    voxel_size_yx = voxel_size[1] # la résolution dans le plan devrait être toujours x/y indep
    cell_mask = cell["cell_mask"]
    nuc_mask = cell["nuc_mask"]
    rna_coord = cell["rna_coord"]
    foci_coord = cell["foci"]
    smfish = cell["smfish"]

    features, features_names = compute_features(cell_mask= cell_mask, nuc_mask= nuc_mask, ndim= 3, rna_coord= rna_coord, smfish= smfish, foci_coord= foci_coord, voxel_size_yx= voxel_size_yx,
        centrosome_coord=None,
        compute_distance=True,
        compute_intranuclear=True,
        compute_protrusion=True,
        compute_dispersion=True,
        compute_topography=True,
        compute_foci=True,
        compute_area=True,
        return_names=True)
    
    header = ["id", "AcquisitionId"] + features_names
    data = np.append([0, acquisition_id], features)
    data = data.reshape([1,-1])
    datashape_ref = DataFrame.newframe_Cell()
    new_Cell = pd.DataFrame(data= data, columns= header)
    dataOp.check_samedatashape(new_Cell, datashape_ref) # Ensure datashape stability along different runs

    return new_Cell



def _get_varname(var):
    #To be used  within a function.
    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var][0]



def get_datetime():
    return dt.datetime.now().strftime("%Y%m%d %H-%M-%S \n")





#TODO : Create output.py
def print_parameters(path_out, *parameters, printDateTime= True):
    """
    Print parameters into a .txt file
    
    Parameters
    ----------
        path_out : str
            full_path to saving location.
        *parameters : args(int, float, str)
            parameters to print into text file.
    """

    check_parameter(path_out = (str))
    if path_out[len(path_out)-1] == '/' : path_out += 'parameters.txt'
    elif path_out[len(path_out)-4 : len(path_out)] != '.txt' : path_out += '.txt'

    parameter_file = open(path_out, "w")

    #Header
    parameter_file.write("PARAMETERS\n")
    parameter_file.write("\n############\n")

    if printDateTime:
        datetime = get_datetime()
        parameter_file.write(datetime)
        parameter_file.write("############\n\n")
    
    lines= []
    for parameter in parameters :
        name = _get_varname(parameter)
        lines += ["{0} : {1}\n".format(name, parameter)]
    parameter_file.writelines(lines)

    parameter_file.close()


def print_dict(dic, path_out):
    check_parameter(path_out = (str), dic = (dict))
    if path_out[len(path_out) -1] == '/' : path_out += 'dic.txt'
    elif path_out[len(path_out)-4 : len(path_out)] != '.txt' : path_out += '.txt'

    lines = []
    for elmt in dic:
        lines += "{0} : {1}\n".format(elmt, dic[elmt])


    dict_file = open(path_out, "w")
    dict_file.writelines(lines)
    dict_file.close()