"""
This sub module contains decorators function to apply on top of plot functions
"""


import pandas as pd
import CustomPandasFramework.PBody_project.update as update

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
