import pandas as pd
import os
import re
import numpy as np
import CustomPandasFramework.DataFrames as DataFrame
import CustomPandasFramework.operations as dataOp
from bigfish.stack import check_parameter,read_image
from bigfish.classification import compute_features

def get_Input(input_path, channels_list) :
    """ Returns a panda dataFrame filled with input folder info. To learn about expected data informations/shape see newFrame_Input.   
     
     Parameters
    ----------
    input_path : string.
        full path to the input folder as a string.
    channels_list : List[str].
        List of strings indicating the different channels of an acquisition. Thoses strings MUST be contained in the filenames otherwise they will be filtered.

    Returns
    -------
    Input : pd.Dataframe
        Dataframe containing informations on files within the input folder. See newFrame_Input for more informations.
        
        """     
     
    #Integrity checks
    check_parameter(input_path = (str), channels_list = (list))
    for string in channels_list : check_parameter(string = (str))
    
    #filename
    Input = newFrame_Input()
    filenames = os.listdir(input_path)
    Input["filename"] = filenames

    #Channel labelling
    for channel in channels_list :
        Input.loc[
        Input.filename.str.match(".*({0}).*".format(channel)), ["channel"]
        ]= channel
    Input = Input.dropna(subset=["channel"])

    #extension
    extensions = []
    ext_regex = "\.\w*$"
    for file in Input.loc[:,"filename"] :
        extensions += re.findall(ext_regex, file)
    Input.loc[:, "file extension"] = extensions

    #root filename
    root_filenames = []
    for line in Input.loc[:,["filename","channel"]].itertuples(index= False) :
        filename, channel = line
        regex = "({0})|(\.\w*$)".format(channel)
        root_filenames += [re.sub(regex, "", filename)]
    Input.loc[:, "root filename"] = root_filenames
    ##Acquisition index
    #Integrity
    channel_num = len(channels_list)
    Input_groupbyroot = Input.sort_values(["filename"]).value_counts(subset= "root filename")
    if not all(Input_groupbyroot == channel_num) : raise Exception("Some acquisition are missing at least one channel. Please check the completeness of files placed in input.")
    
    #Computing acquisition index
    Input_acquisitionindex = Input_groupbyroot.reset_index(drop = False)
    Input_acquisitionindex = Input_acquisitionindex.rename(columns={0 : "channel-complete"})
    Input_acquisitionindex = Input_acquisitionindex.reset_index(drop = False)
    Input_acquisitionindex = Input_acquisitionindex.rename(columns={"index" : "acquisition index"})
    Input_acquisitionindex = Input_acquisitionindex.drop(["channel-complete"], axis= 1)
    Input = Input.drop("acquisition index", axis= 1)
    Input = pd.merge(Input, Input_acquisitionindex,
            how= "left", left_on= "root filename", right_on= "root filename")
    Input = Input.sort_values(["acquisition index","channel"]).reset_index(drop= True)

    return Input



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
    regex = "(\w*)--"
    res = re.findall(regex,root)[0]
    return res

def newFrame_Input() :
    """Returns an empty pandas DataFrame with expected input data shape. 
    This frame is used for navigating through files put in the input folder during the segmentation process and isn't meant to be stored, therefor it does not contains any key.
        
    Returns
    -------
    
        new_Input : pd.Dataframe
            filename : full name of files within the input folder.
            channel : name of the channel such as dapi or egfp.
            root filename : filename without channel and extension, refering to the acquisition. In our project should be like : RNAspecies-wellID--sk1fk1fl1. It does NOT contains channel and extension informations.
            file extension : extension of the file.
            acquisition : Int. Computed, should be one per root filename, it aims at grouping different channels into the same acquisition.
            
    """
    new_Input = pd.DataFrame({
        "filename" : [],
        "channel" : [],
        "root filename" : [],
        "file extension" : [],
        "acquisition index" : []

        })
    
    return(new_Input)




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