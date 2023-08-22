"""
This submodules groups all function related to box plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CustomPandasFramework.PBody_project.update as update
from CustomPandasFramework.integrity import MissingColumnError
from .utils import save_plot
from .decorators import rotate_xAxis_label

@rotate_xAxis_label
def count_pbody_per_Cell(Pbody: pd.DataFrame, xlabel= None, ylabel= "pbodies per cell", title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    if "rna name" not in Pbody.columns : raise MissingColumnError("`rna name` wasn't found in Pbody frame, consider using update.AddRnaName.")

    group = Pbody.groupby(["rna name", "CellId"])["id"].count().groupby("rna name", level= 0, axis= 0).apply(list)

    
    box_plot(data= group, ylabel= ylabel, labels= group.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)


@rotate_xAxis_label
def count_rna_per_Cell(detection_view: pd.DataFrame, xlabel= None, ylabel= "rna spot per cell", title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    data_mean = detection_view.loc["rna",:].groupby("rna name", level= 0, axis=0)["count"].apply(list)
  
    box_plot(data= data_mean, ylabel= ylabel, labels= detection_view.index.get_level_values(1).unique(), xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)

@rotate_xAxis_label
def cell_number(Acquisition: pd.DataFrame, xlabel= None, ylabel= "cell number", title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :

    """
    1 box per gene
    """
    Df_Acquisition = Acquisition.groupby(["rna name"])["cell number"].apply(list).sort_index()

    box_plot(data= Df_Acquisition, labels= Df_Acquisition.index, xlabel= xlabel, ylabel= ylabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)

@rotate_xAxis_label
def count_Malat_per_Cell(CellCellular_cycle, xlabel= None, ylabel= "malat spot per cell", title= None, reset= True, close= False, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """

    df = CellCellular_cycle.copy().dropna(subset= ["count_in_nuc", "count_in_cyto"])
    df["total"] = df["count_in_nuc"] + df["count_in_cyto"]

    data_mean = df.groupby(by="rna name", axis= 0, level= 0)["total"].apply(list)
    
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)


@rotate_xAxis_label
def dapi_signal(Cell: pd.DataFrame, projtype= 'mean', summarize_type= 'mean', integrated_signal = False,
                 xlabel= None, ylabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene.
    Integrated signal --> True : multiply value by nucleus area.
    """

    if 'rna name' not in Cell.columns : raise MissingColumnError("rna name wasn't found in Cell DataFrame. Consider using update.AddRnaName")

    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    
    if integrated_signal:
        Cell["integrated signal"] = Cell[X] * Cell["nucleus area (nm^2)"]
        X = "integrated signal"
    
    if ylabel == None : ylabel = X

    dataframe = Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).mean().loc[:,[X]]
    dataframe.sort_values("rna name").reset_index(drop= False)
    data_mean = Cell.groupby("rna name")[X].apply(list)
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)


## Base plot ##

def raw_data(Cell: pd.DataFrame, Acquisition: pd.DataFrame, column_name, xlabel= None, ylabel= None, title= None, reset= False, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])

    Df_Acquisition = Join_Cell.loc[:,["rna name","AcquisitionId", column_name]].groupby(["rna name","AcquisitionId"]).mean().loc[:,[column_name]]

    Df_Acquisition = Df_Acquisition.sort_values("rna name")
    Df_Acquisition = Df_Acquisition.reset_index(drop= False)
    data_mean = Df_Acquisition.groupby("rna name")[column_name].apply(list)
    
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)


def box_plot(data: np.ndarray, xlabel= None, ylabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    Default plot for box plots.

    Parameters
    ----------
        data : sequence[float]
        
        **kargs :
            color

    """

    if reset : fig = plt.figure(figsize=(20,10))
    else : fig = plt.gcf()

    plt.boxplot(data, **kargs)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output=path_output, ext=ext)
    
    if show : plt.show()
    if close : plt.close()

    return fig