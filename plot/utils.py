import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
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


def gene_bar_plot(rna_list: 'list[str]', values: 'list[float]', errors: 'list[float]',
                 title: str=None, xlabel:str= None, ylabel:str= None,
                path_output= None, ext='png', show= True, close= True, legend: 'list[str]'= None, width = 0.8, error_width = 3) :

    #Exception thrower
    is_listoflist = False
    if show and not close : raise ValueError("If show is True then close shoud not be False as closing the graph windows will close the plot.")

    #multi data set
    if type(values[0]) == list : 
        if type(errors[0]) != list : raise TypeError("When passing several bar sets to plot, it is expected that several errors sets be passed as well")
        is_listoflist = True
    if type(errors[0]) == list : 
        if type(values[0]) != list : raise TypeError("When passing several errors sets to plot, it is expected that several bar sets be passed as well")
        is_listoflist = True

    #len list matches
    if not is_listoflist :
        if not(len(rna_list) == len(values) == len(errors)) : raise ValueError("rna, values and errors lengths must match")
    else :
        #Set lengths match
        if not(len(values) == len(errors)) and legend == None : raise ValueError("value sets and errors sets lengths must match")
        elif not(len(values) == len(errors) == len(legend)) : raise ValueError("value sets and errors sets and legend lengths must match")
        #Data lengths match
        for set in range(0,len(values)) :
            if not(len(rna_list) == len(values[set]) == len(errors[set])) : raise ValueError("values and errors lengths must match")

    
    #Init plot
    color_list = ['red','blue','green','orange','purple','brown','cyan'] * (round(len(rna_list)/7) + 1)
    fig = plt.figure(figsize= (20,10))
    ax = fig.gca()

    #Case when several bars are plotted for each genes
    if is_listoflist :
        bar_number = len(values)
        length = width/bar_number
        error_width /= bar_number
        abs = np.arange(0,len(values[0]))
        barshift = np.arange(-(width/2 - length/2),(width/2), step = length)
        color_list = color_list[:len(barshift)]
        assert len(barshift) == len(values), "barshift : {0} ; values : {1}".format(len(barshift), len(values))
        for bar_set, error_set, shift, color, label in zip(values, errors, barshift, color_list, legend) :
            X = abs - shift
            if legend != None : ax.bar(X, bar_set, yerr= error_set, capsize= error_width, color= color, width= length, align= 'center', label= label)
            else : ax.bar(X, bar_set, yerr= error_set, capsize= error_width, color= color, width= length, align= 'center')
        if legend != None : ax.legend()
    
    #Case one bar per gene
    else :
        plt.bar(rna_list, values, yerr= errors, capsize= error_width, color= color_list[:len(rna_list)], width= width, align= 'center')
    
    plt.axis(ymin= 0)
    
    #Spacing
    plt.xticks(range(len(rna_list)))
    xticks = ax.get_xticks()
    ax.set_xticks(xticks, labels= rna_list, rotation= 90)
    ax.axis(xmin=-0.5, xmax= len(rna_list)+0.5, ymin= 0)
    fig.subplots_adjust(bottom= 2/fig.get_size_inches()[1])
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    if close : plt.close()

    return fig



def gene_violin_plot(gene_list: 'list[str]', values: 'list[float]', errors: 'list[float]'):
    """
    
    """
    pass



def histogram(data: 'list[float]', xlabel= 'distribution', ylabel= 'count', path_output= None, show = True, close= True, ext= 'png', title: str = None, bins= 500, ticks_number= 21, **axis_boundaries) :
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
    for arg in axis_boundaries:
        if arg not in ('xmin','xmax','ymin','ymax') : raise ValueError("Incorrect argument in axis_boundaries, expected arguments are 'xmin','xmax','ymin','ymax' not {0}.".format(arg))
    if type(data) != np.ndarray:
        data = np.array(data)
    
    #Plot
    data = data[~np.logical_or(np.isnan(data),np.isinf(data))] # ~ = not
    fig = plt.figure(figsize= (20,10))
    hist = plt.hist(data, bins= bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

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
        if xmax > 0 : last_tick = ceil(xmax)
        else : last_tick = 0
        if xmin < 0 : first_tick = floor(xmin)
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
    

def format_array_scientific_notation(array) :
    res = map(functools.partial(np.format_float_scientific, precision= 2), array)
    return list(res)
