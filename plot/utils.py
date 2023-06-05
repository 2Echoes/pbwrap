import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter
from math import floor, ceil
import functools

def from_label_get_centeroidscoords(label: np.ndarray):
    """
    Returns
    --------
      centroid : dict{"label": list, "centroid-n": list} 
        n should be replace with 1,2.. according to the axe you wish to access."""

    check_parameter(label = (np.ndarray))
    centroid = regionprops_table(label, properties= ["label","centroid"])
    return centroid


def get_colors_list(size) :
    color_list  = list(mcolors.CSS4_COLORS.keys())
    for color in ['white', 'snow', 'whitesmoke', 'ivory', 'floralwhite', 'ghostwhite', 'seashell', 'linen', 'honeydew', 'aliceblue', 'mintcream'] :
        color_list.remove(color)
    length = len(color_list) - 1
    index = np.linspace(0,length,size).round().astype(int)
    color_list = np.array(color_list)

    return color_list[index]




def set_axis_ticks(axis:tuple, x_ticks_number:int = None, y_ticks_number: int= None) :
    """
    Set 'ticks_number' ticks on the plot, ticks are spaced regularly using min/max value from axis tuple.

    Parameters
    ----------
        axis : tuple (xmin, xmax, ymin, ymax)
        ticks_number : int
    """
    if not isinstance(axis, (tuple, list)) : raise TypeError("axis paremeter should be a tuple or list. It is a {0}".format(type(axis)))
    if len(axis) != 4 : raise ValueError("axis parameter should be a list containing 4 float-like : xmin, xmax, ymin, ymax.")

    xmin,xmax,ymin,ymax = axis

    #X axis
    if x_ticks_number != None :
        if xmax > 0 : last_tick = round(xmax)
        else : last_tick = 0
        if xmin < 0 : first_tick = round(xmin)
        else : first_tick = 0   
        x_ticks = np.linspace(first_tick,last_tick,x_ticks_number)
        if all(np.abs(x_ticks) > 1) : x_ticks = np.round(x_ticks)
        else : x_ticks = np.round(x_ticks, decimals= 3)
        x_ticks[0] = xmin
        x_ticks[x_ticks_number-1] = xmax
        if any(x_ticks >= 10000) : x_label = format_array_scientific_notation(x_ticks)
        elif all(x_ticks < 1) and all(x_ticks > -1) : x_label = format_array_scientific_notation(x_ticks)
        else : x_label = x_ticks
        xlocs, xlabels = plt.xticks(x_ticks,x_label)

    else : xlocs, xlabels = None,None

    #Y axis
    if y_ticks_number != None :
        last_tick = ceil(ymax) 
        y_ticks = y_ticks = np.linspace(0,last_tick,y_ticks_number)
        y_ticks[0] = floor(xmin)
        y_ticks[y_ticks_number-1] = ymax
        if any(y_ticks >= 10000) : x_label = format_array_scientific_notation(y_ticks)
        elif all(y_ticks< 1) and all(y_ticks > -1) : y_label = format_array_scientific_notation(y_ticks)
        else : y_label = y_ticks
        ylocs, ylabels = plt.xticks(y_ticks, y_label)
    else : ylocs, ylabels = None,None

    return(xlocs,xlabels,ylocs,ylabels)

    




def save_plot(path_output, ext):
    """Save the plot.

    Parameters
    ----------
    path_output : str
        Path to save the image (without extension).
    ext : str or List[str]
        Extension used to save the plot. If it is a list of strings, the plot
        will be saved several times.
    
    Code from BigFish package.
    BSD 3-Clause License

    Copyright © 2020, Arthur Imbert
    All rights reserved.
    """
    # add extension at the end of the filename
    if ext == None : ext ='png'
    extension = "." + ext
    if extension not in path_output:
        path_output += extension

    # save the plot
    if isinstance(ext, str):
        # add extension at the end of the filename
        extension = "." + ext
        if extension not in path_output:
            path_output += extension
        plt.savefig(path_output, format=ext)
    elif isinstance(ext, list):
        for ext_ in ext:
            # add extension at the end of the filename
            extension = "." + ext_
            if extension not in path_output:
                path_output += extension
            plt.savefig(path_output, format=ext_)
    else:
        Warning("Plot is not saved because the extension is not valid: "
                "{0}.".format(ext))
        

def round_up(x) :
    x = ((x // 10)+1)*10
    return x


def round_up_bis(x,digit) :
    x = ceil(x/10**digit)*10**digit
    return x


def truncate(x, digit) :
    x = x * 10**digit
    x = ceil(x)
    x = x/10**digit
    return x
    


def format_array_scientific_notation(array: np.ndarray, precision = 2) :
    """
    Format an iterable of float into scientific notation using numpy scientific notation.
    """
    res = map(functools.partial(auto_format_float_scientific, precision= precision), array)
    return list(res)



def auto_format_float_scientific(number: float, precision: int):
    """
    Format each element from an iterable of float with more than 5 digits into scientific notation using numpy scientific notation.
    Never formats 0.

    """
    
    if number == 0 : res = 0
    elif len(str(float)) < 5 :
        res = number
    else : res = np.format_float_scientific(number,precision=precision)
    
    return res

