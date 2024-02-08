"""
This submodules groups all function related to scatter plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools as itools

from ..utils import check_parameter
from .utils import multi_plot_positions
from .utils import save_plot


## Base plot ##

def plot(X: np.ndarray, Y: np.ndarray, xlabel= None, ylabel= None, title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    #TODO does not function bc of label handling.
    """
    Default plot for points plots..query("`rna name` in ['NF1', 'PABPC1']")

    Parameters
    ----------
        data : sequence[float]
        
        **kargs :
            color

    """
    #auto kargs
    if "marker" not in kargs :
        kargs["marker"] = '.'
    if "ls" not in kargs :
        kargs['ls'] = ''


    if reset : fig = plt.figure(figsize=(20,10))
    else : fig = plt.gcf()

    plt.plot(X,Y, **kargs)
    if "label" in kargs :
        plt.legend()


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if path_output != None : save_plot(path_output=path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


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

    if reset : fig = plt.figure(figsize=(20,20))
    else : fig = plt.gcf()


    #Auto kargs :
    if not 'edgecolors' in kargs :
        kargs['edgecolors'] = 'black'
    
    if "label" in kargs and "color" in kargs :
        label,color = kargs["label"], kargs["color"]
        del kargs["label"], kargs["color"]
        set_legend(labels= label, colors= color, **kargs)
        kargs["label"], kargs["color"] = label,color
        del label,color


    plt.scatter(X,Y, **kargs)
    
    plt.axis('tight')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if path_output != None : save_plot(path_output=path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


    return fig

def vertical_scatter(
        ax: plt.Axes, distributions, labels=None, colors=None, xlabel=None, ylabel=None, title=None, y_axis= None,
        showmeans=True, show_std=True, multi_distributions= False
                ):
    
    #Check parameters
    check_parameter(ax = plt.Axes, xlabel = (type(None), str), ylabel = (type(None), str), title = (type(None),str), y_axis= (type(None),tuple,list), multi_distributions= bool)
    if multi_distributions :
        for distrib in distributions :
            check_parameter(distrib = (list, tuple, np.ndarray, pd.DataFrame, pd.Series))

    if type(labels) == type(None) : pass
    elif len(distributions) != len(labels) : 
        print("distributions : ", distributions)
        print("labels : ", labels)
        raise ValueError("Length of labels and distributions must match.")

    if type(colors) == type(None): pass
    elif len(distributions) != len(colors) and not multi_distributions: raise ValueError("Length of colors and distributions must match.")
    
    if type(y_axis) == type(None) : pass
    elif len(y_axis) != 2 : raise ValueError("y_axis is expected to be of length 2 : (ymin, ymax).")
    
    if multi_distributions :
        if type(colors) != type(None) :
            if len(colors) == len(distributions) :
                new_colors = []
                for color,distrib in zip(colors, distributions) :
                    new_colors.extend([color]*len(distrib))
                colors = new_colors

        max_individual_violin_number = max([len(distrib) for distrib in distributions]) + 1 #Is the maximum number of group plotted for one set.
        positions, ticks_positions = multi_plot_positions(distributions)
        distributions = list(itools.chain(*distributions))
        assert len(distributions) == len(positions), "AssertionError : multi_distributions wrongly flattened : positions : {0}, distributions {1}".format(len(positions), len(distributions))

        #colors
        if type(colors) != type(None) :
            if len(colors) == len(distributions) : pass
            else : raise ValueError("Length of colors must either match length of distributions or the number of element in distributions")
    
    else :
        positions = np.arange(1, len(distributions) + 1)
        ticks_positions = np.arange(1, len(distributions) + 1)
        max_individual_violin_number = 1
    
    new_positions = []
    new_colors = []
    for sample_index, group in enumerate(distributions):
        new_positions.extend([positions[sample_index]]* len(group))
        new_colors.extend([colors[sample_index]]* len(group))

    #Plot
    scatter_plot = ax.scatter(
        y = distributions,
        x = new_positions,
        s = 10,
        c = new_colors,
        alpha= 0.8
        )

    if type(labels) == type(None) :
        labels = np.arange(1, len(distributions) + 1)
    xticks = ax.set_xticks(ticks_positions, labels=labels)

    ax.set_xlim(0.25, len(labels) * max_individual_violin_number + 0.75)
    
    if type(y_axis) != type(None) :
        axis = list(ax.axis())
        if y_axis[0] != None : axis[2] = y_axis[0]
        if y_axis[1] != None : axis[3] = y_axis[1]
        ax.axis(axis)
    else : axis = ax.axis()
    
    if showmeans :
        means = [np.mean(distrib) for distrib in distributions]
        ax.scatter(positions, means, c= colors, s= 40, marker= 'D', linewidths=0.5, edgecolors='black', alpha= 0.5)
    
    if show_std :
        std = [np.std(distrib) for distrib in distributions]
        ax.errorbar(positions, means, std, capsize = 3, ecolor= 'black', fmt = 'None', elinewidth= 1)

    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)
    if type(ylabel) != type(None) : ax.set_ylabel(ylabel)
    if type(title) != type(None) : ax.set_title(title)

    return ax

def set_legend(labels, colors, column_number = 3, loc = None, **kargs) :
    df = pd.DataFrame(data = {"label" : labels, "color" : colors})
    df = df.value_counts(subset=["label", "color"]).reset_index(drop= False)
    for label, color in zip(df["label"], df["color"]) :
        plt.scatter([],[],label = label, color = color, **kargs)
    
    plt.legend(ncol= column_number, loc=loc,)