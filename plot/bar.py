"""
This submodules groups all function related to bar plots making from base plot to result plots.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pbwrap.data.getdata as gdata
import CustomPandasFramework.PBody_project.update as update
import CustomPandasFramework.PBody_project.views as views
from .utils import save_plot, check_parameter
from itertools import chain
from matplotlib.container import BarContainer
from matplotlib.patches import Rectangle

def threshold(Acquisition: pd.DataFrame, rna_list:'list[str]' = None, path_output= None, show = True, close= True, ext= 'png', title = None) :
    
    #Computing RNA mean threshold and var :
    if rna_list == None : rna_list =  gdata.from_Acquisition_get_rna(Acquisition)
    elif type(rna_list) == str : rna_list = [rna_list]
    
    threshold_list = [] 
    std_list = []
    for rna in rna_list :
        threshold_list += [Acquisition[Acquisition["rna name"] == rna].loc[:,"RNA spot threshold"].mean()]
        std_list += [Acquisition[Acquisition["rna name"] == rna].loc[:,"RNA spot threshold"].std()]

    fig = gene_bar_plot(rna_list, threshold_list, errors= std_list, title= title, ylabel= "mean threshold", path_output= path_output, ext=ext, show=show, close= close)

## P-Body ##
def P_body_detect_inside_nucleus(Cell: pd.DataFrame, path_output= None, show = True, close = True, ext= 'png', title: str = None) :
    """
    Plot a 2 bar graphs : one bar for number of pbodies detected inside nucleus and one for number of pbodies detected in cytoplasm.
    """

    Number_in_cytoplasm = Cell.loc[:, "count pbody in cytoplasm"].sum()
    std_in_cytoplasm = Cell.loc[:, "count pbody in cytoplasm"].std()
    Number_in_nucleus = Cell.loc[:, "count pbody in nucleus"].sum()
    std_in_nucleus = Cell.loc[:, "count pbody in nucleus"].std()

    ylabel = "Count"

    gene_bar_plot(["P-bodies in cytoplasm", "P-Bodies in nucleus"], [Number_in_cytoplasm, Number_in_nucleus], [std_in_cytoplasm, std_in_nucleus],ylabel= ylabel, path_output=path_output, show=show, close=close, ext=ext, title=title)



## DAPI ##
def CellularCycle_classification(Cell_list, plate_list, cellular_cycle_classifier= update.from_IntegratedSignal_spike_compute_CellularCycleGroup,
                                 title: str='Celular_cycle repartition from plate to plate', xlabel:str= 'plate', ylabel:str= 'proportion', path_output= None, ext='png', show= True, close= True, width = 0.8, error_width = 3) :
    """
    Computes Cellular cycle repartion on a each plates of plate_list using ``cellular_cycle_classifier`` and creates bar graph to compare overall distributions.

    Parameters
    ----------
        Cell_list : list[pd.DataFrames]
        plate_list : list['str']
        cellular_cycle_classifier : func 
            Must add 'cellular_cycle' key to a Cell DataFrame.
    """
    
    Cell_frames = (cellular_cycle_classifier(Cell) for Cell in Cell_list)
    view_list = (views.CellularCycle_distribution_view(frame) for frame in Cell_frames)
    dictionary_list = [views.CellularCycle_distribution_overallmeasure(view) for view in view_list]
    
    #Unpacking results
    g1_values = [di['g1_mean'] for di in dictionary_list]
    g1_std = [di['g1_std'] for di in dictionary_list]
    g2_values = [di['g2_mean'] for di in dictionary_list]
    g2_std = [di['g2_std'] for di in dictionary_list]
    S_values = [di['s_mean'] for di in dictionary_list]
    S_std = [di['s_std'] for di in dictionary_list]

    gene_bar_plot(rna_list= plate_list, values= [g1_values, g2_values, S_values], errors= [g1_std, g2_std, S_std], legend= ['g1', 'g2', 's'],
                  title=title, xlabel=xlabel, ylabel= ylabel, path_output=path_output, ext=ext, show=show, close=close, width=width, error_width=error_width)

    



def DapiSignal_InfValue(Acquisition:pd.DataFrame, Cell:pd.DataFrame, max_value: float, gene_list=None, projtype= 'mean', summarize_type = 'median', path_output= None, show = True,close= True, ext= 'png', title: str = None):
    """
    byGenes_barplot
    projtype : "MIP" or "MEAN"
    summarize_type : "median" or "mean"

    Standard deviation is calculated from an Acquisition point of view

    """

    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    std_data= []
    mean_data = []
    join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
    join_frame = join_frame.drop(axis= 0, index= join_frame[join_frame["pbody number"] == 0].index)
    if gene_list == None : gene_list = gdata.from_Acquisition_get_rna(Acquisition)
    for gene in gene_list : 
        gene_Cell = join_frame [join_frame["rna name"] == gene]
        cell_proportion_under_value = np.array([])

        for acquisition in gene_Cell.value_counts(subset= "AcquisitionId").index :
            acquisition_Cell = gene_Cell.query("AcquisitionId == {0}".format(acquisition))
            cell_under_value = (len(acquisition_Cell.query("{0} <= {1}".format(X,max_value))))
            total_cell = (len(acquisition_Cell))
            cell_proportion_under_value = np.append(cell_proportion_under_value, (cell_under_value / total_cell))

        mean_data.append(cell_proportion_under_value.mean())
        std_data.append(cell_proportion_under_value.std())


    gene_bar_plot(gene_list, mean_data, std_data, title=title, ylabel=X, path_output=path_output, ext=ext, show=show, close= close)


## Spots plots ##

def spots_per_cell(Acquisition: pd.DataFrame, Cell: pd.DataFrame, spot_type = 'rna', rna_list: 'list[str]' = None, path_output= None, show = True, close=True, ext= 'png', title = "RNA per cell") :
    
    #Determining spot type
    if spot_type.upper() == 'RNA' : column = "rna number"
    elif spot_type.upper() == 'MALAT1' : 
        column = "malat1 number"
        Cell["malat1 number"] = Cell.loc[:, "malat1 spots in cytoplasm"] + Cell.loc[:, "malat1 spots in nucleus"]
    else : raise ValueError("spot type shoud either be 'rna' or 'malat1'. {0} was given".format(spot_type))


    if rna_list == None : rna_list =  gdata.from_Acquisition_get_rna(Acquisition)
    elif type(rna_list) == str : rna_list = [rna_list]

    Cell = gdata.from_rna_get_Cells(rna= rna_list, Cell= Cell, Acquisition= Acquisition)
    Cell["rna number"] = Cell["nb_rna_out_nuc"] + Cell["nb_rna_in_nuc"]

    #Computing mean values and std
    threshold_list = []
    std_list = []
    for rna in rna_list :
        threshold_list += [Cell[Cell["rna name"] == rna].loc[:, column].mean()]
        std_list += [Cell[Cell["rna name"] == rna].loc[:, column].std()]

    #plot
    fig = gene_bar_plot(rna_list, threshold_list, std_list, title= title, path_output= path_output, ext=ext, show=show, close= close)


def cytoRNA_proportion_in_pbody(detection_view:pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, close=True, ext= 'png', title = "Nucleus RNA proportion inside P-bodies"):

    if gene_list == None : gene_list = detection_view.index.get_level_values(1).unique()
    std_list= []
    mean_list = []
    for gene in gene_list :
        group = detection_view.loc[("rna",gene), ["count_in_cyto", "count_in_Pbody"]]
        hasPbodyidx = group.query("count_in_Pbody == 0").index
        group = group.drop(hasPbodyidx, axis= 0)
        group["proportion"] = group["count_in_Pbody"]/ group["count_in_cyto"]
        mean = group["proportion"].mean()
        std = group["proportion"].std()
        mean_list.append(mean)
        std_list.append(std)

    fig = gene_bar_plot(gene_list, mean_list, std_list, title= title,ylabel="Proportion of cytoplasmic rna detected inside P-bodies.",
                         path_output= path_output, ext=ext, show=show, close= close)



def RNA_proportion_in_pbody(detection_view: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, close= True, ext= 'png', title = "RNA proportion inside P-bodies"):

    if gene_list == None : gene_list = detection_view.index.get_level_values(1).unique()
    
    std_list= []
    mean_list = []
    for gene in gene_list : 
        group = detection_view.loc[("rna",gene), ["count", "count_in_Pbody"]]
        
        group["proportion"] = group["count"]/ group["count_in_Pbody"]
        hasPbodyidx = group.query("count_in_Pbody == 0").index
        group = group.drop(hasPbodyidx, axis= 0)
        mean = group["proportion"].mean()
        std = group["proportion"].std()
        mean_list.append(mean)
        std_list.append(std)

    fig = gene_bar_plot(gene_list, mean_list, std_list, title= title,ylabel="Proportion of cytoplasmic rna detected inside P-bodies.",
                         path_output= path_output, ext=ext, show=show, close= close)



def RNApercentage_in_nucleus(detection_view: pd.DataFrame, gene_list: 'list[str]' = None, plot_in_and_out_bars= True, ylabel= 'rna proportion inside nucleus', path_output= None, show = True, close= True, ext= 'png', title = None) :

    if gene_list == None : gene_list = detection_view.index.get_level_values(1).unique()
    std_list= []
    mean_list = []
    for gene in gene_list : 
        group = detection_view.loc[("rna",gene), ["count_in_nuc", "count"]]
        group["proportion"] = group["count_in_nuc"]/ group["count"]
        mean = group["proportion"].mean()
        std = group["proportion"].std()
        mean_list.append(mean)
        std_list.append(std)

    fig = gene_bar_plot(gene_list, mean_list, std_list, title= title, ylabel=ylabel,
                         path_output= path_output, ext=ext, show=show, close= close)

def total_cell_number(Acquisition: pd.DataFrame,xlabel=None, ylabel= "Cell number", path_output= None, show = True, close= True, ext= 'png', title = "Computed cell number") :
    """
    1 bar per gene
    """

    Df = Acquisition.groupby(["rna name"])["cell number"].sum()
    gene_bar_plot(Df.index, Df, title= title, ylabel= ylabel, xlabel=xlabel, show=show, close=close, path_output=path_output, ext=ext )


#######################
###### BASE PLOT ######
#######################
    
def bar_plot(ax: plt.Axes, data, errors=None,
            labels=None, colors=None, xlabel=None, ylabel=None, title=None, y_axis=None,
            multi_bar_plot= False, 
            **kwargs
            ) :
    
    """
    Parameters
    ----------

    data : list[bar_height] (or list[list[bar_height]] for multi_bar_plot = True)
    errors : list[float-like] : one per bar.
    color : list[color] 1 per bar number or 1 per group of bars.
    multi_bar_plot : pass True if data is listof list.
    """
    
    #Parameters check
    check_parameter(ax = plt.Axes, 
                    data = (list, tuple, np.ndarray, pd.DataFrame, pd.Series, pd.Index), errors = (list, tuple, np.ndarray, pd.DataFrame, pd.Series, pd.Index, type(None)), labels = (list, tuple, np.ndarray, pd.DataFrame, pd.Series, pd.Index, type(None)), colors= (list, tuple, np.ndarray, pd.DataFrame, pd.Series, pd.Index, type(None)),
                    xlabel= (str,type(None)), ylabel= (str,type(None)), title= (str,type(None)),
                    y_axis= (tuple, type(None)),
                    multi_bar_plot= bool)
    
    #Catch KWARGS
    #logscale
    if kwargs.get('logscale') :
        logscale=True
        del kwargs['logscale']
    else: logscale= False


    if multi_bar_plot :
        for iterable in data : check_parameter(iterable = (list, tuple, np.ndarray, pd.DataFrame, pd.Series, pd.Index))
        if type(errors) != type(None) : 
            for iterable in errors : check_parameter(iterable = (list, tuple, np.ndarray, pd.DataFrame, pd.Series, pd.Index))

    number_of_set = len(data)

    if type(labels) != type(None) :
        if len(labels) != number_of_set : raise ValueError("length of labels must match length of data.")
      

    #Layout
    if multi_bar_plot :

        if type(colors) != type(None) :
            if len(colors) == len(data) :
                new_colors = []
                for color,distrib in zip(colors, data) :
                    new_colors.extend([color]*len(distrib))
                colors = new_colors

        max_bar_number_per_set = max(len(set) for set in data)# + 1 # +1 for spacing bewteen groups of bars
        positions, xticks_positions = multi_plot_positions(data)
        data = list(chain(*data))
        if type(errors) != type(None) : errors = list(chain(*errors))
        assert len(data) == len(positions), "AssertionError : multi_bars_data wrongly flattened : positions : {0}, distributions {1}".format(len(positions), len(data))
    
    else : 
        max_bar_number_per_set = 1
        positions = np.arange(1, len(data) + 1)
        xticks_positions = np.arange(1, len(data) + 1)
    
    ##### From there data is a flat iterable

    #Plot
        kwargs.setdefault('capsize', 10)
    bar_plot = ax.bar(
        x= positions, 
        height= data,
        yerr= errors,
        width= 0.75,
        align= 'center',
        color= colors,
        label= 'data_bars',
        **kwargs
           )
    
    # bar_handles(bar_plot, colors)

    #x axis labels
    if type(labels) == type(None) : labels = np.arange(1,len(xticks_positions) + 1)
    xticks = ax.set_xticks(xticks_positions, labels=labels)
    ax.set_xlim(0.25, max(positions) + 0.75)
    
    if logscale :
        ax.set_yscale('log')

    #colors   
    if type(colors) != type(None) :
        if len(colors) != len(data) and len(colors) != number_of_set : raise ValueError("Length of colors must either match length of data or the number of element in data")
    
    if type(y_axis) != type(None) :
        axis = list(ax.axis())
        if y_axis[0] != None : axis[2] = y_axis[0]
        if y_axis[1] != None : axis[3] = y_axis[1]
        ax.axis(axis)
    else : axis = ax.axis()

    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)
    if type(ylabel) != type(None) : ax.set_ylabel(ylabel)
    if type(title) != type(None) : ax.set_title(title)

    return ax
    


def gene_bar_plot(rna_list: 'list[str]', values: 'list[float]', errors: 'list[float]'=None,
                 title: str=None, xlabel:str= None, ylabel:str= None,
                path_output= None, ext='png', show= True, close= True, legend: 'list[str]'= None, width = 0.8, error_width = 3) :
    
    """
    Base Plot for bar graph, rna_list can also take a list of float to display on X axis.
    """

    #Exception thrower
    is_listoflist = False
    if show and not close : raise ValueError("If show is True then close shoud not be False as closing the graph windows will close the plot.")

    #multi data set
    if type(values[0]) == list : 
        if type(errors[0]) != list : raise TypeError("When passing several bar sets to plot, it is expected that several errors sets be passed as well")
        is_listoflist = True
    if type(errors) == type(None) : pass
    elif type(errors[0]) == list : 
        if type(values[0]) != list : raise TypeError("When passing several errors sets to plot, it is expected that several bar sets be passed as well")
        is_listoflist = True

    #len list matches
    if not is_listoflist :
        if type(errors) == type(None) :
            if not(len(rna_list) == len(values)) : raise ValueError("rna and values lengths must match")

        else:
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


def multi_plot_positions(distributions) :

    max_individual_violin_number = max([len(distrib) for distrib in distributions]) + 1#Is the maximum number of violin plotted for one set.

    positions = []
    ticks_positions = []
    for distrib_number, distrib in enumerate(distributions) :
        positions.extend(list(
            np.arange(1, len(distrib) + 1) + (distrib_number * max_individual_violin_number) if len(distrib) > 1 
            else [distrib_number * max_individual_violin_number + (max_individual_violin_number-1)/2 + 1]
        ))

        ticks_positions.append(
            distrib_number * max_individual_violin_number + (len(distrib)-1)/2 + 1 if len(distrib) > 1
            else distrib_number * max_individual_violin_number + (max_individual_violin_number-1)/2 + 1
        )

    return positions, ticks_positions

def _make_bar_handles(color_list) :
    return [Rectangle(xy= (0,0), width= 0.75, height = 1, color= color) for color in color_list]