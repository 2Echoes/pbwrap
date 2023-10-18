import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import CustomPandasFramework.PBody_project.get as get
from .utils import get_colors_list, get_markers_generator, hide_overlapping_annotations, save_plot


## Base plot ##

def _Layout_Quantif_plots(gene_list, gene_outlier_dict:dict =None, xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    Makes a plot containing up to 3 G1G2 subplots.

    Returns
    -------
        fig, ax, modified kargs.

    """
    if len(gene_list) == 0 : raise ValueError("gene_list is empty.")

    if 'figsize' in kargs :
        figsize = kargs['figsize']
    else :
        figsize = (45, 15)

    #Main Plot settings
    if reset : fig, ax = plt.subplots(nrows= 1, ncols= 3, figsize= figsize)
    else : fig = plt.gcf()
    if type(title) != type(None) : plt.suptitle(title, fontsize= 20, fontweight= 'bold')
    if type(xlabel) != type(None) : fig.text(0.5, 0.04, xlabel, ha='center', fontsize= 30)

    #Color set
    if "color" in kargs :
        color_list = kargs["color"]
        if len(color_list) == 1 : color_list *= len(gene_list)
        elif isinstance(color_list, str) : color_list = [color_list]*len(gene_list)
    else : 
        color_list = get_colors_list(len(gene_list))

    #linewidth set
    if "linewidths" in kargs :
        linewidths_list = kargs["linewidths"]
        if len(linewidths_list) == 1 : linewidths_list *= len(gene_list)
        elif isinstance(linewidths_list, (float,int)) : linewidths_list = [linewidths_list]*len(gene_list)
    else : linewidths_list = [1]*len(gene_list)

    #edgecolor set
    if "edgecolors" in kargs :
        edgecolors_list = kargs["edgecolors"]
        if len(edgecolors_list) == 1 : edgecolors_list *= len(gene_list)
        elif isinstance(edgecolors_list, str) : edgecolors_list = [edgecolors_list]*len(gene_list)
    else : edgecolors_list = ['black']*len(gene_list)

    if type(gene_outlier_dict) != None :

        for key in ['Df', 'number', 'pk'] :
            if key not in gene_outlier_dict.keys() : raise KeyError('{0} key not found in gene_outlier_dict')
        Df = gene_outlier_dict['Df'].sort_values(['rna name'])
        number = gene_outlier_dict['number']
        pk = gene_outlier_dict['pk']

        red_genes = get.get_GenesWithlessThanXCells(Df.query("cellular_cycle == 'g1'"),number,pk)
        red_genes += get.get_GenesWithlessThanXCells(Df.query("cellular_cycle == 'g2'"),number,pk)
        red_index = [gene_list.index(gene) for gene in red_genes]
        for idx in red_index :
            edgecolors_list[idx] = 'red'
            linewidths_list[idx] = 1.2

    kargs["color"] = color_list
    kargs["linewidths"] = linewidths_list
    kargs["edgecolors"] = edgecolors_list
    kargs["alpha"] = 0.7

    return fig, ax, kargs

def _G1G2_main_legend_layout(fig) :
    #Main plot legend
    handle1 = plt.scatter([],[], color= "white", linewidths= 1.5, edgecolor='black', label= 'Genes with more \nthan 100 cells computed')
    handle2 = plt.scatter([],[], color= "white", linewidths= 2, edgecolor='red', label= 'Genes with less \nthan 100 cells computed')
    handles = [handle1, handle2]
    labels = ['Genes with more than 100 cells computed', 'Genes with less than 100 cells computed']
    fig.legend(handles, labels, loc='upper left', prop={'size': 10})

    return fig

def G1G2_plot(Data : pd.Series, plot_X_equal_Y_line= True, cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
              title= None, legend= True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs):

    """
    Makes a G1 versus G2 plot. Expected Data Input is a pandas Series with multi-index ['rna name', 'cellular_cycle', 'CellId'].

    Output a scatter plot with a point per gene with x value = Serie cell average for 'cellular_cycle' = 'g1' and y value = Serie cell average for 'cellular_cycle' = 'g2'.

    Parameters
    ----------
        Data : pd.Series
            Expected Data Input is a pandas Series with multi-index ['rna name', 'cellular_cycle', 'CellId'].
        ...

        kargs :
            edgecolors (list)
            linewidths (list)
            color (list)
            axis [xmin, xmax, ymin, ymax] or ['square', 'tight', ...]
    """ 
    
    if not isinstance(Data, pd.Series) : raise TypeError("Data argument must be a pd.Series type.")
    level = Data.index.nlevels
    if level != 2 and level != 3 : raise IndexError("Data argument index is expected to be a multi-index with nlevels = 2 or 3 (CellId optional): ['gene name', 'cellular_cycle', ]. {0} levels were found".format(Data.index.nlevels))
    for index_name in ['rna name', 'cellular_cycle'] : 
        if index_name not in Data.index.names : raise KeyError("{0} was not found in Data index.".format(index_name))
    
    gene_list = Data.sort_index().index.get_level_values(0).unique()

    if reset : plt.figure(figsize=(20,20))

    kargs_copy = kargs.copy()
    if 'color' in kargs :
        colors = iter(kargs['color'])
        del kargs_copy['color']
    else : colors = iter(get_colors_list(len(gene_list)))
    
    if not 'edgecolors' in kargs :
        kargs['edgecolors'] = ['black'] * len(gene_list)
    else : del kargs_copy["edgecolors"]
    if not 'linewidths' in kargs :
        kargs['linewidths'] = [1] * len(gene_list)
    else : del kargs_copy['linewidths']

    markers_gen = get_markers_generator()
    annotation_list = []
    for gene, lw, edgecolor in zip(gene_list, kargs["linewidths"], kargs["edgecolors"]) :
        marker = next(markers_gen)
        color =  next(colors)
        DF = Data.loc[gene,:]
        index_lvl0 = DF.index.get_level_values(0).unique()
        if cellular_cycle_x in index_lvl0  and level == 3 : g1_mean = DF.loc[cellular_cycle_x,:].mean()
        elif cellular_cycle_x in index_lvl0  and level == 2 : g1_mean = DF.loc[cellular_cycle_x].mean()
        else : g1_mean = 0
        if cellular_cycle_y in index_lvl0 and level == 3 : g2_mean = DF.loc[cellular_cycle_y,:].mean()
        elif cellular_cycle_y in index_lvl0  and level == 2 : g2_mean = DF.loc[cellular_cycle_y].mean()
        else : g2_mean = 0
        plt.scatter(x= g1_mean, y= g2_mean, color = color, label= gene, linewidths=lw, marker=marker, edgecolors= edgecolor, s= 60, **kargs_copy)
        annotation_list.append(plt.text(x= g1_mean*0.98, y = g2_mean*1.01, s= gene, size= 10))
    
    if legend : plt.legend(ncols= 4)
    hide_overlapping_annotations(*annotation_list)

    if 'axis' in kargs :
        plt.axis(kargs['axis'])
    else :
        xmin,xmax,ymin,ymax = plt.axis('square')
        plt.axis([0,xmax,0,ymax])
    if plot_X_equal_Y_line :
        X = [0, xmax]
        plt.plot(X,X,'b')

    plt.xlabel(cellular_cycle_x.upper())
    plt.ylabel(cellular_cycle_y.upper())
    if type(title) != type(None) : plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()