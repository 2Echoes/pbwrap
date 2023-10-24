import pandas as pd
import datetime as dt
import re, inspect
import CustomPandasFramework.operations as dataOp
from bigfish.stack import check_parameter,read_image




def get_images_as_gen(path_input: str, Input: pd.DataFrame, acquisition_list: 'list[int]', channels_list: 'list[str]'= None, z_min= None, z_max= None) :
    """ Open images from an acquisition within the Input DataFrame using bigfish open method. 
        Outputs as a generator of images ordered by rootfilename then by channel_name (A->Z) .
    
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

    #Integrity checks)
    check_parameter(path_input = (str), Input = (pd.DataFrame), acquisition_list = (list, int), channels_list = (list, str, type(None)))

    if type(channels_list) == str : channels_list = [channels_list]
    if type(acquisition_list) == int : acquisition_list = [acquisition_list]
    if channels_list == None :
         channels_list = Input.values_count(subset= "channel").index.tolist()


    # images = []
    for acquisition in acquisition_list :
        index = Input.query("`acquisition index` == {0} and channel in {1}".format(acquisition, channels_list)).sort_values(["filename", "channel"]).index
        for fileindex in index:
            filename = Input.at[fileindex, "filename"]
            images = read_image(path= path_input + filename)
            if images.ndim == 3 and (z_min != None or z_max != None):
                if z_max == None : z_max = images.shape[0]
                if z_min == None : z_min = 0
                images = images[z_min: z_max,:,:]
            yield images



def get_images_as_list(path_input: str, Input: pd.DataFrame, acquisition_list: 'list[int]', channels_list: 'list[str]'= None, z_min= None, z_max = None) :
    """ Open images from an acquisition within the Input DataFrame using bigfish open method. 
        Outputs as a generator of images ordered by rootfilename then by channel_name (A->Z) .
    
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

    #Integrity checks)
    check_parameter(path_input = (str), Input = (pd.DataFrame), acquisition_list = (list, int), channels_list = (list, str, type(None)))

    if type(channels_list) == str : channels_list = [channels_list]
    if type(acquisition_list) == int : acquisition_list = [int]
    if channels_list == None :
         channels_list = Input.values_count(subset= "channel").index.tolist()


    # images = []
    images = []
    for acquisition in acquisition_list :
        index = Input.query("`acquisition index` == {0} and channel in {1}".format(acquisition, channels_list)).sort_values(["filename", "channel"]).index
        for fileindex in index:
            filename = Input.at[fileindex, "filename"]
            image = read_image(path= path_input + filename)
            if image.ndim == 3 and (z_min != None or z_max != None):
                if z_max == None : z_max = image.shape[0]
                if z_min == None : z_min = 0
                image = image[z_min: z_max,:,:]

            images += [image]

    return images

def get_acquisition_num(Input: pd.DataFrame) :
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

    acquisition_num = len(Input.value_counts(subset=["acquisition index"]))
    return acquisition_num


def get_rootfilename(acquisition_index, Input_frame: pd.DataFrame):
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

    idx = Input_frame.query("`acquisition index` == {0}".format(acquisition_index)).index
    res = Input_frame.at[idx[0], "root filename"]
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

def from_rootfilename_get_Acquisitionid(Acquisition: pd.DataFrame, rootfilename:str):
    return Acquisition.query("rootfilename == '{0}'".format(rootfilename)).reset_index(drop=True).at[0,"id"]


def from_Acquisition_get_rna(Acquisition: pd.DataFrame) :
    return list(Acquisition.value_counts(subset= "rna name").sort_index().index)



def from_rna_get_Cells(rna: 'list[str]', Cell: pd.DataFrame, Acquisition: pd.DataFrame) -> pd.DataFrame :
    """
    Returns sub-table from Cell containing only cells which rna name (from Acquisition) matches one in 'rna'.
    Also 'rna name' column is added to returned Cell sub-table.
    """


    if type(rna) == str : rna = [rna]

    join_frame = dataOp.keep_columns(Dataframe= pd.merge(Cell, Acquisition, how= 'left', left_on= 'AcquisitionId', right_on= 'id'),
                                     columns= ["rna name"] + list(Cell.columns))
    drop_index = join_frame.query('`rna name` not in {0}'.format(rna)).index
    join_frame = join_frame.drop(axis= 0, index= drop_index)

    return join_frame



def _get_varname(var):
    #To be used  within a function.
    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var][0]



def get_datetime():
    return dt.datetime.now().strftime("%Y%m%d_%H-%M-%S")