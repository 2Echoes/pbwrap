"""
This sub module contains decorators function to apply on top of plot functions
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import CustomPandasFramework.PBody_project.update as update
from .utils import save_plot

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