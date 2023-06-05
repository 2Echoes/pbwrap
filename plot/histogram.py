"""
This submodules groups all function related to histogram making from base plot to result plots.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CustomPandasFramework.PBody_project.update as update
from .utils import set_axis_ticks, save_plot
from .decorators import use_g1_g2_grouping

########################
##### RESULT PLOTS #####
########################

def dapi_signal(Cell, projtype= 'MIP', summarize_type = 'median', path_output= None, show = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    From the Maximum Intensity Projection (MIP) or Mean projection computes the histogram of integrated signal within cells (signal*nucleus area).
    Signal can be chosen from mean or median.

    Parameters
    ----------

        projtype : "MIP" or "MEAN"
        summarize_type : "median" or "mean"
        axis_boundaries : kwargs
            boundaries for x and y axes. Expected None or at least one the following ('xmin'=x, 'xmax'=X, ymin='y',ymax='Y')

    """
    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    dapi_signal = Cell.loc[:,X] * Cell.loc[:,"nuc_area"]
    histogram(dapi_signal, xlabel="Dapi signal (Nucleus area * {0})".format(X), ylabel= "count", path_output=path_output, show=show, ext=ext, title=title, bins=bins, **axis_boundaries)


def dapi_density(Cell, projtype= 'MIP', summarize_type = 'median', path_output= None, show = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    From the Maximum Intensity Projection (MIP) or Mean projection computes the histogram of dapi density within cells (signal/nucleus area).
    Signal can be chosen from mean or median.

    Parameters
    ----------

        projtype : "MIP" or "MEAN"
        summarize_type : "median" or "mean"
        axis_boundaries : kwargs
            boundaries for x and y axes. Expected None or at least one the following ('xmin'=x, 'xmax'=X, ymin='y',ymax='Y')

    """
    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    dapi_signal = Cell.loc[:,X] / Cell.loc[:,"nuc_area"]
    histogram(dapi_signal, xlabel="Dapi signal ({0}/ Nucleus area)".format(X), ylabel= "count", path_output=path_output, show=show, ext=ext, title=title, bins=bins, **axis_boundaries)



def malat_count(Cell, location= 'nucleus', path_output= None, show = True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    Histogram of malat spots detected in nucleus, cytoplasm or both.
    location : 'nucleus', 'cytoplasm' or 'cell' (cell = nuc + cytoplasm)
    projtype : "MIP" or "MEAN"
    """
    Cell = update.from_nucleus_malat_proportion_compute_CellullarCycleGroup(Cell, 0.5)
    if location.upper() == "NUCLEUS" : X = "malat1 spots in nucleus"
    elif location.upper() == "CYTOPLASM" or location.upper() == "CELL" : X = "malat1 spots in cytoplasm"
    else : raise ValueError("Incorrect value for location parameter. Should be one of the following : 'nucleus', 'cytoplasm' or 'cell'.")
    dapi_signal = Cell.loc[:, [X, "Cellular_cycle (malat proportion)"]]
    g1 = dapi_signal.query("`Cellular_cycle (malat proportion)` == 'g1'").index
    g2 = dapi_signal.query("`Cellular_cycle (malat proportion)` == 'g2'").index
    if location.upper() == "CELL" : dapi_signal += Cell.loc[:, "malat1 spots in nucleus"]
    
    histogram(dapi_signal.loc[:, X], color= 'green',xlabel=X, ylabel= "count", show=show, path_output=path_output, ext=ext, close=close, title=title, bins=bins, **axis_boundaries)


def in_nuc_malat_proportion(Cell, path_output= None, show = True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    Histogram of malat proportion detected inside nucleus.
    """
    proportion = Cell.loc[:, "malat1 spots in nucleus"] / (Cell.loc[:, "malat1 spots in nucleus"] + Cell.loc[:, "malat1 spots in cytoplasm"])
    histogram(proportion, xlabel="malat spots proportion in nucleus", ylabel= "count", path_output=path_output, show=show, close=close, ext=ext, title=title, bins=bins, **axis_boundaries)





def RawData(DataFrame:pd.DataFrame, variable_name:str, color='blue', label:str=None, path_output= None, show = True, reset= True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries):
    """Basic hist plot for graph requiring the distribution just as it appears in CellDataFrame"""

    data = DataFrame.loc[:,variable_name]
    histogram(data, xlabel= variable_name, color=color, label=label, ylabel= "count", path_output=path_output, show=show, reset=reset, close= close, ext=ext, title=title, bins=bins, **axis_boundaries)






#######################
###### BASE PLOT ######
#######################


def histogram(data: 'list[float]', color= 'blue', label:str = None, xlabel= 'distribution', ylabel= 'count', path_output= None, show = True, close= True, reset= True, ext= 'png', title: str = None, bins= 500, ticks_number= 21, disable_axis= False, **axis_boundaries) :
    """Base function for histograms plotting.
    
    Parameters
    ----------
        data : list
            data distribution
        xlabel : str
            label to plot on x axis
        ylabel : str
            label to plot on y axis
        title : str
            title to plot on top of graph
        axis_boundaries : boundaries for x and y axes. Expected None or at least one the following ('xmin'=x, 'xmax'=X, ymin='y',ymax='Y')
    """
    #Value errors
    if not close and show :
        raise ValueError("if you want to keep the graph open you must not show() during the call of this function")
    if type(data) != np.ndarray:
        data = np.array(data)
    
    #Plot
    data = data[~np.logical_or(np.isnan(data),np.isinf(data))] # ~ = not
    if reset : fig = plt.figure(figsize= (20,10))
    else : fig = plt.gcf()
    
    hist = plt.hist(data, bins= bins, color= color, label= label, alpha= 0.5)
    if not disable_axis : 
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if label != None : plt.legend()

        #Axis boundaries
        ax = fig.gca()
        axis_bound = hist_set_axis_boundaries(ax, data, hist, **axis_boundaries)
        #Axis ticks
        set_axis_ticks(axis_bound, ticks_number)


        if title != None : plt.title(title)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    if close : plt.close()


def hist_set_axis_boundaries(ax, data=None, hist=None, **axis_boundaries) :
    """
    Auto or manual set of histogram boundaries.

    Parameters
    ----------
        ax : figure ax can be get from figure.gca() -> get curret axis
        axis_boundaries : boundaries for x and y axes. Expected None or at least one the following ('xmin'=x, 'xmax'=X, ymin='y',ymax='Y')
        Data and Hist parameters should be given for auto set. (hist is get from hist = plt.hist(data))

    Returns
    -------
        axis = [xmin,xmax,ymin,ymax]

    """

    if 'xmin' in axis_boundaries : xmin = axis_boundaries['xmin']
    else : xmin = 0
    if 'xmax' in axis_boundaries : xmax = axis_boundaries['xmax']
    else : xmax = data.max()
    if 'ymin' in axis_boundaries : ymin = axis_boundaries['ymin']
    else : ymin = 0
    if 'ymax' in axis_boundaries : ymax = axis_boundaries['ymax']
    else : ymax = np.array(hist[0]).max()
    axis = [xmin,xmax,ymin,ymax]
    ax.axis(axis)
    return axis
