"""
This submodules groups all function related to histogram making from base plot to result plots.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CustomPandasFramework.PBody_project.update as update
from .utils import set_axis_ticks, save_plot, get_colors_list
from math import sqrt, floor



########################
##### RESULT PLOTS #####
########################


def dapi_signal(Cell, projtype= 'mean', summarize_type = 'mean', path_output= None, show = True, ext= 'png', title: str = None, bins= 500, auto_bins= True, alpha=0.5, **axis_boundaries) :
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

    dapi_signal = Cell.loc[:,X] * Cell.loc[:,"nucleus area (nm^2)"]
    title += "   {0} cells".format(len(dapi_signal))
    if len(dapi_signal) < 100 and auto_bins :bins = 20


    if "xlabel" not in axis_boundaries :
        xlabel = "Dapi signal (Nucleus area * {0})".format(X)
    else : 
        xlabel = axis_boundaries["xlabel"]
        del axis_boundaries["xlabel"]

    if "ylabel" not in axis_boundaries :
        ylabel = "count"
    else : 
        ylabel = axis_boundaries["ylabel"]
        del axis_boundaries["ylabel"]


     

    distribution = histogram(dapi_signal, xlabel=xlabel, ylabel= ylabel, path_output=path_output, show=show, ext=ext, title=title, bins=bins,alpha=alpha, **axis_boundaries)
    return distribution

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



def malat_count(CellCycle_view, location= 'nucleus', path_output= None, show = True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    Histogram of malat spots detected in nucleus, cytoplasm or both.
    location : 'nucleus', 'cytoplasm' or 'cell' (cell = nuc + cytoplasm)
    projtype : "MIP" or "MEAN"
    """


    df = CellCycle_view.copy().dropna(subset= ["count_in_nuc", "count_in_cyto"])
    if location.upper() == "NUCLEUS" : X = "count_in_nuc"
    elif location.upper() == "CYTOPLASM" or location.upper() == "CELL" : X = "count_in_cyto"
    else : raise ValueError("Incorrect value for location parameter. Should be one of the following : 'nucleus', 'cytoplasm' or 'cell'.")
    dapi_signal = df.loc[:, X]
    if location.upper() == "CELL" : dapi_signal += df.loc[:, "count_in_nuc"]
    
    histogram(dapi_signal, color= 'green',xlabel=X, ylabel= "count", show=show, path_output=path_output, ext=ext, close=close, title=title, bins=bins, **axis_boundaries)


def in_nuc_malat_proportion(CellCycle_view, path_output= None, show = True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    Histogram of malat proportion detected inside nucleus.
    """
    df = CellCycle_view.copy().dropna(subset=["count_in_nuc", "count_in_cyto"])
    proportion = df.loc[:, "count_in_nuc"] / (df.loc[:, "count_in_nuc"] + df.loc[:, "count_in_cyto"])
    histogram(proportion, xlabel="malat spots in nucleus proportion", ylabel= "count", path_output=path_output, show=show, close=close, ext=ext, title=title, bins=bins, **axis_boundaries)


def rna_in_pbodies(Pbody: pd.DataFrame, bins=20, auto_bins= True, path_output= None, show = True, reset=True, close = True, ext= 'png', title: str = "Rna spots in P-bodies Genes Distributions") :
    """
    Tile plot of rna spots distribution in P-bodies.
    """

    df = Pbody.copy().set_index(["rna name", "id"]).fillna(0)

    gene_list = df.index.get_level_values(0).unique()
    #tile init
    plot_number = len(gene_list)
    root = sqrt(plot_number)
    if root - floor(root) == 0 :
        n_lin = n_col = root
    elif root - floor(root) < 0.5 :
        n_lin = floor(root)
        n_col = floor(root) + 1
    else :
        n_lin = floor(root) + 1
        n_col = floor(root) + 1


    if reset : fig, ax = plt.subplots(nrows= int(n_lin), ncols= int(n_col), figsize= (int(12*n_col), int(12*n_lin)))
    else : fig, ax = plt.gcf(), plt.gca()

    plot_idx=1
    colors = iter(get_colors_list(plot_number))
    
    #Subplot
    for gene in gene_list :
        plt.subplot(n_lin, n_col, plot_idx)

        plot_idx +=1
        color = next(colors)
        distribution = df.loc[gene, "rna 0nm count"]
        pbody_number = len(distribution)
        
        if auto_bins :
            if pbody_number < 100 : bins = 20
            else : bins = 50

        histogram(distribution, color, bins= bins, y_label = "count", title= "{0} - {1} P-bodies.".format(gene, pbody_number), show= False, close= False, reset= False, edgecolor = 'black', linewidth= 3)

    fig.text(0.5, 0.04, 'Number of rna spots detected in P-bodies', ha='center', fontsize= 65)
    if type(title) != type(None) : plt.suptitle(title, fontsize= 80, fontweight= 'bold')
    if type(path_output) != type(None) : plt.savefig(path_output + "Rna spots in Pbodies tiles plot")
    if show : plt.show()
    if close : plt.close()



def RawData(DataFrame:pd.DataFrame, variable_name:str, color='blue', label:str=None, path_output= None, show = True, reset= True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries):
    """Basic hist plot for graph requiring the distribution just as it appears in CellDataFrame"""

    data = DataFrame.loc[:,variable_name]
    histogram(data, xlabel= variable_name, color=color, label=label, ylabel= "count", path_output=path_output, show=show, reset=reset, close= close, ext=ext, title=title, bins=bins, **axis_boundaries)






#######################
###### BASE PLOT ######
#######################


def histogram(data: 'list[float]', color= 'blue', label:str = None, xlabel= 'distribution', ylabel= 'count', path_output= None, show = True, close= True, reset= True, ext= 'png', title: str = None, bins= 500, alpha=0.5, ticks_number= 21, disable_axis= False, **axis_boundaries) :
    """
    Base function for histograms plotting. Returns data array used for distribution histogram plot.
    
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
    hist = plt.hist(data, bins= bins, color= color, label= label, edgecolor= 'white', lw=1, alpha= alpha)
    if not disable_axis : 
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if label != None : plt.legend()

        if title != None : plt.title(title)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    if close : plt.close()

    return data


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


