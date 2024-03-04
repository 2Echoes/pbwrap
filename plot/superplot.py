import pandas as pd
import matplotlib.pyplot as plt
from .utils import make_color_frame, get_markers_generator

def distribution_super_plot(data, ax, title=None, xlabel=None, ylabel=None) :
    """
    
    Super plot is returned as plt.Axes
    """
    if type(data) != pd.Series : raise TypeError('data argument passed is of type {0} should be "pd.Series".'.format(type(data)))
    if type(ax) != plt.Axes : raise TypeError('ax argument passed is of type {0} should be "plt.Axes".'.format(type(ax)))

    level = len(data.index.names)

    if level == 1 :
        ax = _distribution_lvl1(data, ax)
    
    elif level == 2 :
        ax = _distribution_lvl2(data, ax)

    elif level == 3 :
        ax = _distribution_lvl3(data, ax)
    
    if type(title) != type(None) : ax.set_title(title)
    if type(xlabel) != type(None) : ax.set_title(title)
    if type(ylabel) != type(None) : ax.set_title(title)


class LevelError(IndexError) :
    pass

def _distribution_lvl1(data: pd.Series, ax: plt.Axes) :

    multi_index = data.index.names
    if len(multi_index) != 1 : raise LevelError("_distribution_lvl1 was called but multi-index dimension does not match.")


def _distribution_lvl2(data: pd.Series, ax : plt.Axes) :
    
    multi_index = data.index.names
    if len(multi_index) != 2 : raise LevelError("_distribution_lvl2 was called but multi-index dimension does not match.")
    data = data.sort_index()
    measure = data.name

    index_distribution_lvl2 = data.index.get_level_values(1).unique()
    colors = make_color_frame(labels= index_distribution_lvl2)



def _distribution_lvl3(data: pd.Series, ax: plt.Axes) :
    """

    yaxis = measure_value

    level1 = x-axis (axis label)
    level2 = distribution colors (legend)
    level3 = smaller scatter points within distribution whith different shapes (legend)

    """
    multi_index = data.index.names
    if len(multi_index) != 3 : raise LevelError("_distribution_lvl3 was called but multi-index dimension does not match.")
    data = data.sort_index()
    measure = data.name

    #Levevel 1&2
    data_distribution_lvl2 = data.reset_index(drop=False).groupby(multi_index[:2])[measure].apply(list)
    ax = _distribution_lvl2(data_distribution_lvl2, ax)

    #Level 3
    index_lvl3 = data.index.get_level_values(2).unique()
    marker_list = list(get_markers_generator(index_lvl3))
    marker_frame = pd.Series(data= marker_list, index= index_lvl3)

    #TODO Récupérer les couleurs utilisées au lvl 2, ainsi que les positions. Class ??