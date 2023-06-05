"""
This submodules groups all function related to scatter plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.update as update
import matplotlib.pyplot as plt
from ..quantification.CurveAnalysis import simple_linear_regression
from .utils import save_plot, get_colors_list



def Malat_inNuc_asDapiIntensity(Cell: pd.DataFrame, projtype = 'MIP', summarize_type= 'median', out = False, plot_linear_regression= False,
                                path_output= None, show = True, ext= 'png', title = None) :
    """
    
    Scatter plot computed as Signal from DAPI channel (X-axis) VS Malat spots count in nucleus.
    Signal can be chosen from mean or median and computed from the Maximum Intensity Projection (MIP) or Mean projection.
    
    """



    #Select projection type from Cell Data
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype shoud either be 'mip' or 'mean', it is {0}.".format(projtype))

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    #Select malat spots in nucleus or in cytoplasm
    if out : Y = 'malat1 spots out nucleus'
    else : Y = 'malat1 spots in nucleus'

    #Fetching data
    Cell["SignalArea"] = Cell[X] * Cell["nuc_area"]
    data = Cell.loc[:, ["SignalArea", Y]].sort_values("SignalArea")

    #Plot
    fig = plt.figure(figsize=(20,10))
    X_values, Y_values = np.array(data.iloc[:,0]), np.array(data.iloc[:,1])
    plt.plot(X_values, Y_values, 'r.', label= "Experimental Data")

    if plot_linear_regression:
        slope, intercept = simple_linear_regression(X_values,Y_values)
        regression = slope* X_values + intercept
        plt.plot(X_values, regression, 'b', label= "Linear regression \n{0}x + {1}".format(slope,intercept))
        plt.legend()

    ax = fig.gca()
    ax.axis(xmin= 0, ymin= 0)
    if title != None : plt.title(title)
    plt.xlabel(X)
    plt.ylabel(Y)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()



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
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
    
    color_list = ( )






## Base plot ## 
def scatter_plot(X: np.ndarray, Y: np.ndarray, xlabel= None, ylabel= None, title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
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


    plt.scatter(X,Y, **kargs)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if show : plt.show()
    if close : plt.close()
    if path_output != None : save_plot(path_output=path_output, ext=ext)

    return fig