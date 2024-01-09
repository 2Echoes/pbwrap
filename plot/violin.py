"""
This submodules groups all function related to violin plots making from base plot to result plots.
"""
import numpy as np
import matplotlib.pyplot as plt
from pbwrap.utils import check_parameter


def violin_plot(
        ax: plt.Axes, distributions, labels=None, colors=None, xlabel=None, ylabel=None, title=None, y_axis= None,
        linewith = 1, line_color = 'black', 
        vertical_plot=True, showmedians= True, showextrema = True
                ) :

    """
    Basis for violin plots.
    """

    #Check parameters

    check_parameter(ax = plt.Axes, xlabel = (type(None), str), ylabel = (type(None), str), vertical_plot= bool, title = (type(None),str), linewith = int, line_color= str, y_axis= (type(None),tuple,list))
    
    if type(labels) == type(None) : pass
    elif len(distributions) != len(labels) : raise ValueError("Length of labels and distributions must match.")

    if type(colors) == type(None) : pass
    elif len(distributions) != len(colors) : raise ValueError("Length of colors and distributions must match.")
    
    if type(y_axis) == type(None) : pass
    elif len(y_axis) != 2 : raise ValueError("y_axis is expected to be of length 2 : (ymin, ymax).")
    #Plot

    violin_plot = ax.violinplot(
        distributions, 
        vert= vertical_plot, 
        showmedians=showmedians,
        showextrema= showextrema
        )

    if type(labels) == type(None) :
        labels = np.arange(1, len(distributions) + 1)
    
    xticks = ax.set_xticks(np.arange(1, len(distributions) + 1), labels=labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    if type(colors) != type(None) :
        for violin, color in zip(violin_plot['bodies'], colors) :
            violin.set_facecolor(color)
            violin.set_alpha(0.6)

    for collection_name in ['cbars', 'cmins', 'cmaxes', 'cmedians'] :
        collection = violin_plot.get(collection_name)
        if type(collection) != type(None) :    
            collection.set_color(line_color)
            collection.set_linewidth(linewith)
    
    if type(y_axis) != type(None) :
        axis = list(ax.axis())
        if y_axis[0] != None : axis[2] = y_axis[0]
        if y_axis[1] != None : axis[3] = y_axis[1]
        ax.axis(axis)

    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)
    if type(ylabel) != type(None) : ax.set_ylabel(ylabel)
    if type(title) != type(None) : ax.set_title(title)

    return violin_plot