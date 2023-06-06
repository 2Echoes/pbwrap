"""
This submodules groups all function related to box plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CustomPandasFramework.PBody_project.update as update
from .utils import save_plot


def count_pbody_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= None, ylabel= "pbodies per cell", title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])

    Df_Acquisition = Join_Cell.loc[:,["rna name","AcquisitionId", "pbody number"]].groupby(["rna name","AcquisitionId"]).mean().loc[:,["pbody number"]]

    Df_Acquisition = Df_Acquisition.sort_values("rna name")
    Df_Acquisition = Df_Acquisition.reset_index(drop= False)
    data_mean = Df_Acquisition.groupby("rna name")["pbody number"].apply(list)
    
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)



def count_rna_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= None, ylabel= "rna spot per cell", title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    Join_Cell["count_rna_per_cell"] = (Join_Cell["nb_rna_out_nuc"] + Join_Cell["nb_rna_in_nuc"])

    Df_Acquisition = Join_Cell.loc[:,["rna name","AcquisitionId", "count_rna_per_cell"]].groupby(["rna name","AcquisitionId"]).mean().loc[:,["count_rna_per_cell"]]

    Df_Acquisition = Df_Acquisition.sort_values("rna name")
    Df_Acquisition = Df_Acquisition.reset_index(drop= False)
    data_mean = Df_Acquisition.groupby("rna name")["count_rna_per_cell"].apply(list)
    
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)


def count_Malat_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= None, ylabel= "malat spot per cell", title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    Join_Cell["count_malat_per_cell"] = (Join_Cell["malat1 spots in nucleus"] + Join_Cell["malat1 spots in cytoplasm"])

    Df_Acquisition = pd.merge(left= Join_Cell.loc[:,["rna name","AcquisitionId", "count_malat_per_cell"]].groupby(["rna name","AcquisitionId"]).mean()["count_malat_per_cell"], 
                              right=Join_Cell.loc[:,["rna name","AcquisitionId", "count_malat_per_cell"]].groupby(["rna name","AcquisitionId"]).std()["count_malat_per_cell"]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={'count_malat_per_cell_x' : 'mean', 'count_malat_per_cell_y' : 'std'})

    Df_Acquisition = Df_Acquisition.sort_values("rna name")
    Df_Acquisition = Df_Acquisition.reset_index(drop= False)
    data_mean = Df_Acquisition.groupby("rna name")["mean"].apply(list)
    
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)



def dapi_signal(Cell: pd.DataFrame, Acquisition: pd.DataFrame, projtype= 'mean', summarize_type= 'mean', integrated_signal = False,
                 xlabel= None, ylabel= None, title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene.
    Integrated signal --> True : multiply value by nucleus area.
    """

    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    
    if integrated_signal:
        Join_Cell["integrated signal"] = (Join_Cell[X] * Join_Cell["nucleus area (nm^2)"])
        X = "integrated signal"
    
    if ylabel == None : ylabel = X


    dataframe = Join_Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).mean().loc[:,[X]]
    dataframe.sort_values("rna name").reset_index(drop= False)
    data_mean = dataframe.groupby("rna name")[X].apply(list)
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)


## Base plot ##

def raw_data(Cell: pd.DataFrame, Acquisition: pd.DataFrame, column_name, xlabel= None, ylabel= None, title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])

    Df_Acquisition = Join_Cell.loc[:,["rna name","AcquisitionId", column_name]].groupby(["rna name","AcquisitionId"]).mean().loc[:,[column_name]]

    Df_Acquisition = Df_Acquisition.sort_values("rna name")
    Df_Acquisition = Df_Acquisition.reset_index(drop= False)
    data_mean = Df_Acquisition.groupby("rna name")[column_name].apply(list)
    
    box_plot(data= data_mean, ylabel= ylabel, labels= data_mean.index, xlabel=xlabel, title= title, reset= reset, close=close, show= show, path_output=path_output, ext=ext, **kargs)



def box_plot(data: np.ndarray, xlabel= None, ylabel= None, title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
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
    
    if show : plt.show()
    if close : plt.close()
    if path_output != None : save_plot(path_output=path_output, ext=ext)

    return fig