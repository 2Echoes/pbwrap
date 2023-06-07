"""
This submodules groups all function related to scatter plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.update as update
import matplotlib.pyplot as plt
from ..quantification.CurveAnalysis import simple_linear_regression
from .utils import save_plot, get_colors_list, identity
from .decorators import plot_curve



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


@plot_curve(identity, label = 'y = x', ls= '-', color= 'black')
def count_Malat_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= "Mean", ylabel= "Standard deviation", title= "Malat spots detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
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
    gene_frame = Join_Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
    color_list = pd.DataFrame(columns = ["color"], data = get_colors_list(gene_number))
    gene_frame = pd.concat([gene_frame, color_list], axis= 1).drop(0, axis= 1)
    Df_Acquisition = pd.merge(left= Df_Acquisition, right= gene_frame, how= 'left', left_on= "rna name", right_on= "rna name")

    scatter(X= Df_Acquisition["mean"], Y = Df_Acquisition["std"], xlabel= xlabel, ylabel= ylabel, color= Df_Acquisition["color"], label= list(Df_Acquisition["rna name"]), title=title, reset=reset, close=close, show=show, path_output=path_output, ext=ext, **kargs)


@plot_curve(identity, label = 'y = x', ls= '-', color= 'black')
def count_rna_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= "Mean", ylabel= "Standard deviation", title= "Rna spots detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    Join_Cell["count_rna_per_cell"] = (Join_Cell["nb_rna_out_nuc"] + Join_Cell["nb_rna_in_nuc"])

    Df_Acquisition = pd.merge(left= Join_Cell.loc[:,["rna name","AcquisitionId", "count_rna_per_cell"]].groupby(["rna name","AcquisitionId"]).mean()["count_rna_per_cell"], 
                              right=Join_Cell.loc[:,["rna name","AcquisitionId", "count_rna_per_cell"]].groupby(["rna name","AcquisitionId"]).std()["count_rna_per_cell"]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={'count_rna_per_cell_x' : 'mean', 'count_rna_per_cell_y' : 'std'})
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_frame = Join_Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
    color_list = pd.DataFrame(columns = ["color"], data = get_colors_list(gene_number))
    gene_frame = pd.concat([gene_frame, color_list], axis= 1).drop(0, axis= 1)
    Df_Acquisition = pd.merge(left= Df_Acquisition, right= gene_frame, how= 'left', left_on= "rna name", right_on= "rna name")

    scatter(X= Df_Acquisition["mean"], Y = Df_Acquisition["std"], xlabel= xlabel, ylabel= ylabel, color= Df_Acquisition["color"], label= list(Df_Acquisition["rna name"]), title=title, reset=reset, close=close, show=show, path_output=path_output, ext=ext, **kargs)


def count_pbody_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= "Mean", ylabel= "Standard deviation", title= "P-bodies detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])

    Df_Acquisition = pd.merge(left= Join_Cell.loc[:,["rna name","AcquisitionId", "pbody number"]].groupby(["rna name","AcquisitionId"]).mean()["pbody number"], 
                              right=Join_Cell.loc[:,["rna name","AcquisitionId", "pbody number"]].groupby(["rna name","AcquisitionId"]).std()["pbody number"]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={'pbody number_x' : 'mean', 'pbody number_y' : 'std'})
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_frame = Join_Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
    color_list = pd.DataFrame(columns = ["color"], data = get_colors_list(gene_number))
    gene_frame = pd.concat([gene_frame, color_list], axis= 1).drop(0, axis= 1)
    Df_Acquisition = pd.merge(left= Df_Acquisition, right= gene_frame, how= 'left', left_on= "rna name", right_on= "rna name")

    scatter(X= Df_Acquisition["mean"], Y = Df_Acquisition["std"], xlabel= xlabel, ylabel= ylabel, color= Df_Acquisition["color"], label= list(Df_Acquisition["rna name"]), title=title, reset=reset, close=close, show=show, path_output=path_output, ext=ext, **kargs)




def dapi_signal(Cell: pd.DataFrame, Acquisition: pd.DataFrame, projtype= 'mean', summarize_type= 'mean', integrated_signal = False,
                xlabel= "Mean", ylabel= "Standard deviation", title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    Integrated signal -> True : signal value is multiplied with nucleus area (nm)
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
        Join_Cell["Integrated signal ({0})".format(X)] = (Join_Cell[X] * Join_Cell["nucleus area (nm^2)"])
        X = "Integrated signal ({0})".format(X)

    if title == None : title = X

    Df_Acquisition = pd.merge(left= Join_Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).mean()[X], 
                              right=Join_Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).std()[X]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={ X+'_x' : 'mean', X+'_y' : 'std'})
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_frame = Join_Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
    color_list = pd.DataFrame(columns = ["color"], data = get_colors_list(gene_number))
    gene_frame = pd.concat([gene_frame, color_list], axis= 1).drop(0, axis= 1)
    Df_Acquisition = pd.merge(left= Df_Acquisition, right= gene_frame, how= 'left', left_on= "rna name", right_on= "rna name")

    scatter(X= Df_Acquisition["mean"], Y = Df_Acquisition["std"], xlabel= xlabel, ylabel= ylabel, color= Df_Acquisition["color"], label= list(Df_Acquisition["rna name"]), title=title, reset=reset, close=close, show=show, path_output=path_output, ext=ext, **kargs)



## Base plot ##

def plot(X: np.ndarray, Y: np.ndarray, xlabel= None, ylabel= None, title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    #TODO does not function bc of label handling.
    """
    Default plot for points plots.

    Parameters
    ----------
        data : sequence[float]
        
        **kargs :
            color

    """

    if hasattr(X[0], '__iter__') and hasattr(Y[0],'__iter__') :
        if len(X) != len(Y) : raise ValueError("X and Y iterables must have the same length")
        is_list = True
    elif not hasattr(X[0], '__iter__') and not hasattr(Y[0],'__iter__') :
        is_list = False
    else : raise TypeError("Either both X-elements and Y-elements should be iterables or none")


    if reset : fig = plt.figure(figsize=(20,10))
    else : fig = plt.gcf()

    if is_list :
        if "color" in kargs :
            color = kargs["color"]
            if len(color) != len(X) : raise ValueError("length of color array must match length of data-set array")
            del kargs["color"]

        else : color = get_colors_list(len(X))

        if "label" in kargs :
            label = kargs["label"]
            if len(label) != len(X) : raise ValueError("length of label array must match length of data-set array")
            del kargs["label"]
        else : label = [None] * len(X)
        
        
        if "ls" in kargs :
            ls = kargs["label"]
            if len(ls) != len(X) : raise ValueError("length of LineStyle (ls) array must match length of data-set array")
            del kargs["ls"]
        else : ls = [None] * len(X)

        for set_number, (x,y) in enumerate(zip(X,Y)) :
            plt.plot(x,y, color= color[set_number], ls= ls[set_number], label= label[set_number], **kargs)
        plt.legend()


    else :
        plt.plot(X,Y, **kargs)
        if "label" in kargs :
            plt.legend()


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if show : plt.show()
    if close : plt.close()
    if path_output != None : save_plot(path_output=path_output, ext=ext)

    return fig



def scatter(X: np.ndarray, Y: np.ndarray, xlabel= None, ylabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    Default plot for scatter plots.

    Parameters
    ----------
        data : sequence[float]
        
        **kargs :
            color

    """

    if reset : fig = plt.figure(figsize=(20,10))
    else : fig = plt.gcf()

    if "label" in kargs and "color" in kargs :
        set_legend(kargs["label"], kargs["color"])
        del kargs["label"]

    plt.scatter(X,Y, **kargs)
    plt.axis('tight')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if path_output != None : save_plot(path_output=path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


    return fig


def set_legend(labels, colors, column_number = 3, loc = 'upper left') :
    df = pd.DataFrame(data = {"label" : labels, "color" : colors})
    df = df.value_counts(subset=["label", "color"]).reset_index(drop= False)
    for label, color in zip(df["label"], df["color"]) :
        plt.scatter([],[],label = label, color = color)
    

    plt.legend(ncol= column_number, loc=loc,)