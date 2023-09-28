"""
This sub module contains decorators function to apply on top of plot functions
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import CustomPandasFramework.PBody_project.update as update
from .utils import save_plot, hist_maximum

def use_g1_g2_grouping(plot_function) :
    """
    Add a color grouping for cell in g1(green) and cell in g2(red)

    Parameters
    ----------
        plot_function : Must have Cell, color and reset as argument.
    """

    def inner(Cell: pd.DataFrame, **kargs) :
        if "Cellular_cycle (malat proportion)" not in Cell.columns :
            Cell = update.from_nucleus_malat_proportion_compute_CellullarCycleGroup(Cell)

        g1_index = Cell.query("`Cellular_cycle (malat proportion)` == 'g1'").index
        g2_index = Cell.query("`Cellular_cycle (malat proportion)` == 'g2'").index

        g1kargs = kargs.copy()
        g1kargs["reset"] = True
        g1kargs["close"] = False
        g1kargs["show"] = False
        g1kargs["color"] = 'green'
        g1kargs["alpha"] = 0.5
        g1kargs["disable_axis"] = True
        g1kargs["label"] = 'g1'
        g1kargs["path_output"] = None
        
        kargs["color"] = 'red'
        kargs["alpha"] = 0.5
        kargs["reset"] = False
        kargs["label"] = 'g2'

        plot_function(Cell.loc[g1_index,:],**g1kargs)
        plot_function(Cell.loc[g2_index,:], **kargs)
    return inner


def rotate_xAxis_label(plot_function) :
    """
    Rotate x axis ticks.
    """
    def inner(labels_list: list=None, rotation=90, **kargs) :

        #Plot function and keep it open
        new_kargs = kargs.copy()
        new_kargs["close"] = False
        new_kargs["show"] = False
        new_kargs["path_output"] = None
        plot_function(**new_kargs)

        #xticks rotation
        ax = plt.gca()
        if labels_list != None :
            xticks = ax.get_xticks() 
            plt.xticks(range(len(labels_list)))
        else :
            xticks, labels_list = plt.xticks()
        
        ax.set_xticks(xticks, labels= labels_list, rotation= rotation)

        #back to called parameters
        if "path_output" in kargs :
            if type(kargs["path_output"]) != type(None) : 
                if "ext" not in kargs : kargs["ext"] = 'png'
                save_plot(path_output= kargs["path_output"], ext= kargs["ext"])
        if "show" in kargs :
            if kargs["show"] : plt.show()
        if "close" in kargs :
            if kargs["close"] : plt.close()

    return inner



def plot_curve(math_func, points_number = 100, legend_col_number=3, **curve_kargs) :
    """
    Add a curve to the plot.
    
    Parameters
    ----------
        math_func : function
            A function defining the curve you wish to add. In other words from a np.ndarray X as argument should return a np.ndarray Y such as : Y = math_func(X).
        points_number : int
            number of points plotted in the graph
        **curve_kargs : any
            **kargs used in matplotlib.pyplot.plot()
    
    """
    def decorator (plot_function) :
        def inner(*args, **kargs) :

            new_kargs = kargs.copy()
            new_kargs["close"] = False
            new_kargs["show"] = False
            new_kargs["path_output"] = None

            plot_function(*args, **new_kargs)

            #adding curve to plot
            xmin, xmax, ymin, ymax = plt.axis()
            X = np.linspace(xmin,xmax, points_number)
            Y = math_func(X)
            plt.plot(X,Y,**curve_kargs)
            if "label" in curve_kargs : 
                plt.legend(ncol= legend_col_number)

            #back to called parameters
            if "path_output" in kargs : 
                if "ext" not in kargs : kargs["ext"] = 'png'
                save_plot(path_output= kargs["path_output"], ext= kargs["ext"])
            if "show" in kargs :
                if kargs["show"] : plt.show()
            else : plt.show()
            if "close" in kargs :
                if kargs["close"] : plt.close()
            else: plt.close()

        return inner
    return decorator


def plot_distribution_median(color= 'orange', anotate= False, last_decorator= True, label= 'auto'):
    """
    Distribution function has to return data it used to plot the distribution.

    Parameters
    ----------
        distribution : np.ndarray, numpy compatible iterable
    """
    def decorator(plot_function):
        def inner(*args, **kargs) :
            new_kargs = kargs.copy()
            new_kargs["close"] = False
            new_kargs["show"] = False
            new_kargs["path_output"] = None

            distribution = plot_function(*args, **new_kargs)
            median = np.median(distribution)
            if median > 10000 : median_anot = np.format_float_scientific(median, 3)
            else : median_anot = median 
            plt.axis('tight')
            xmin, xmax, ymin, ymax = plt.axis()
            X = [median] * len(distribution)
            Y = [ymin] *( len(distribution) - 1 )+ [ymax]
            if anotate : plt.annotate("Median value : {0}".format(median_anot),[median + (xmax-xmin)*0.02, ymax]) # xmin + ,
            if label == 'auto' : legend = "Median value : {0}".format(median_anot)
            else : legend = label
            plt.plot(X,Y, color= color, label= legend)

            if last_decorator :
                #back to called parameters
                plt.legend()
                if "path_output" in kargs : 
                    if "ext" not in kargs : kargs["ext"] = 'png'
                    save_plot(path_output= kargs["path_output"], ext= kargs["ext"])
                if "show" in kargs :
                    if kargs["show"] : plt.show()
                else : plt.show()
                if "close" in kargs :
                    if kargs["close"] : plt.close()
                else: plt.close()
            return distribution
        return inner
    return decorator


def plot_distribution_percentile(percentile: float, color= 'red', anotate= False, last_decorator= True, label= 'auto'):
    """
    Distribution function has to return data it used to plot the distribution.

    Parameters
    ----------
        distribution : np.ndarray, numpy compatible iterable
    """
    def decorator(plot_function):
        def inner(*args, **kargs) :
            new_kargs = kargs.copy()
            new_kargs["close"] = False
            new_kargs["show"] = False
            new_kargs["path_output"] = None

            distribution = plot_function(*args, **new_kargs)
            percentile_value = np.percentile(distribution, percentile)
            if percentile_value > 10000 : percentile_value_anot = np.format_float_scientific(percentile_value, 3)
            else : percentile_value_anot = percentile_value
            # plt.axis('tight')
            xmin, xmax, ymin, ymax = plt.axis()
            X = [percentile_value] * len(distribution)
            Y = [ymin] *( len(distribution) - 1 )+ [ymax]
            if anotate : plt.annotate("{1} percentile value : {0}".format(percentile_value_anot, percentile),[percentile_value + (xmax-xmin)*0.02, ymax]) # xmin + ,
            if label == 'auto' : legend= "{1} percentile value : {0}".format(percentile_value_anot, percentile)
            else : legend = label
            
            plt.plot(X,Y, color= color, label= legend)
            if last_decorator :
                plt.legend()
                #back to called parameters
                if "path_output" in kargs : 
                    if "ext" not in kargs : kargs["ext"] = 'png'
                    save_plot(path_output= kargs["path_output"], ext= kargs["ext"])
                if "show" in kargs :
                    if kargs["show"] : plt.show()
                else : plt.show()
                if "close" in kargs :
                    if kargs["close"] : plt.close()
                else: plt.close()
            return distribution
        return inner
    return decorator


def plot_hist_max(color= 'green', anotate= False, last_decorator= True, multiplicator= 1, shift = 0, auto_shift= False, label= 'auto'):
    """
    Distribution function has to return data it used to plot the distribution.

    Parameters
    ----------
        distribution : np.ndarray, numpy compatible iterable
    """
    
    def decorator(plot_function):
        def inner(*args, **kargs) :
            new_kargs = kargs.copy()
            new_kargs["close"] = False
            new_kargs["show"] = False
            new_kargs["path_output"] = None
            if "bins" in new_kargs.keys() : bins = new_kargs["bins"]
            else : bins = new_kargs["bins"]= 50
            distribution = plot_function(*args, **new_kargs)
            if len(distribution) < 100 : bins = 20

            if auto_shift :
                maximum_value = (hist_maximum(np.histogram(distribution, bins=bins)) - np.percentile(distribution, 0.5)/2 )* multiplicator
            else : maximum_value = (hist_maximum(np.histogram(distribution, bins=bins)) + shift )* multiplicator

            
            if maximum_value > 10000 : maximum_value_anot = np.format_float_scientific(maximum_value, 3)
            else : maximum_value_anot = maximum_value
            plt.axis('tight')
            xmin, xmax, ymin, ymax = plt.axis()
            
            X = [maximum_value] * len(distribution)
            Y = [ymin] *( len(distribution) - 1 )+ [ymax]
            if label == 'auto' : legend = "Histogram maximum value : {0}".format(maximum_value_anot)
            else : legend = label
            if anotate : plt.annotate("Histogram maximum value : {0}".format(maximum_value_anot),[maximum_value + (xmax-xmin)*0.02, ymax]) # xmin + ,
            plt.plot(X,Y, color= color, label= legend)
            if last_decorator :
                plt.legend()
                #back to called parameters
                if "path_output" in kargs : 
                    if "ext" not in kargs : kargs["ext"] = 'png'
                    save_plot(path_output= kargs["path_output"], ext= kargs["ext"])
                if "show" in kargs :
                    if kargs["show"] : plt.show()
                else : plt.show()
                if "close" in kargs :
                    if kargs["close"] : plt.close()
                else: plt.close()
            return distribution
        return inner
    return decorator