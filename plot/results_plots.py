import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
import CustomPandasFramework.PBody_project.update as update
from ..quantification.CurveAnalysis import simple_linear_regression
from .utils import save_plot, gene_bar_plot, histogram


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


def hist_RawData(DataFrame:pd.DataFrame, variable_name:str, color='blue', barlabel:str=None, path_output= None, show = True, reset= True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries):
    """Basic hist plot for graph requiring the distribution just as it appears in CellDataFrame"""

    data = DataFrame.loc[:,variable_name]
    histogram(data, xlabel= variable_name, color=color, barlabel=barlabel, ylabel= "count", path_output=path_output, show=show, reset=reset, close= close, ext=ext, title=title, bins=bins, **axis_boundaries)


#########################
## RNA Detection plots ##
#########################

#############
# Bar plots #
#############

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


def cluster_per_cell() :
    pass #TODO

def RNA_in_pbody(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, close= True, ext= 'png', title = "Count of RNA in P-bodies"):

    join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
    join_frame = join_frame.drop(axis= 0, index= join_frame[join_frame["pbody number"] == 0].index)
    if gene_list == None : gene_list = gdata.from_Acquisition_get_rna(Acquisition)
    
    std_list= []
    mean_rna_per_pbody_list = []
    for gene in gene_list : 
        gene_Cell = join_frame [join_frame["rna name"] == gene]
        mean_rna_per_pbody_list += [(gene_Cell.loc[:, "rna spots in pbody"] / gene_Cell.loc[:, "pbody number"]).mean()]
        std_list += [(gene_Cell.loc[:, "rna spots in pbody"] / gene_Cell.loc[:, "pbody number"]).std()]

    #plot
    fig = gene_bar_plot(gene_list, mean_rna_per_pbody_list, std_list, title= title, xlabel="count", path_output= path_output, ext=ext, show=show, close= close)



def cytoRNA_proportion_in_pbody(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, close=True, ext= 'png', title = "Nucleus RNA proportion inside P-bodies"):

    join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
    join_frame = join_frame.drop(axis= 0, index= join_frame[join_frame["pbody number"] == 0].index)
    if gene_list == None : gene_list = gdata.from_Acquisition_get_rna(Acquisition)
    
    std_list= []
    mean_rna_per_pbody_list = []
    for gene in gene_list : 
        gene_Cell = join_frame [join_frame["rna name"] == gene]
        mean_rna_in_pbody = gene_Cell.loc[:,"rna spots in pbody"].mean()
        mean_cyto_rna_number = gene_Cell.loc[:,"nb_rna_out_nuc"].mean()
        mean_rna_per_pbody_list += [mean_rna_in_pbody/mean_cyto_rna_number]        
        std_rna_in_pbody = gene_Cell.loc[:,"rna spots in pbody"].std()
        std_cyto_rna_number = gene_Cell.loc[:,"nb_rna_out_nuc"].std()
        std_list += [std_rna_in_pbody/std_cyto_rna_number] # Comment calculer la std ?


    #plot
    fig = gene_bar_plot(gene_list, mean_rna_per_pbody_list, std_list, title= title,ylabel="cytoplasmic RNA proportion detected inside p-bodies",
                         path_output= path_output, ext=ext, show=show, close= close)



def RNA_proportion_in_pbody(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, close= True, ext= 'png', title = "RNA proportion inside P-bodies"):

    join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
    join_frame = join_frame.drop(axis= 0, index= join_frame[join_frame["pbody number"] == 0].index)
    if gene_list == None : gene_list = gdata.from_Acquisition_get_rna(Acquisition)
    
    # std_list= []
    mean_rna_per_pbody_list = []
    for gene in gene_list : 
        gene_Cell = join_frame [join_frame["rna name"] == gene]

        mean_rna_per_pbody_list += [gene_Cell.loc[:, "rna spots in pbody"].sum() / (gene_Cell.loc[:, "nb_rna_out_nuc"] + gene_Cell.loc[:, "nb_rna_in_nuc"]).sum()]
        # mean_rna_per_pbody_list += [(gene_Cell.loc[:, "rna spots in pbody"] / (gene_Cell.loc[:, "nb_rna_out_nuc"] + gene_Cell.loc[:, "nb_rna_in_nuc"])).mean()]
        # std_list += [(gene_Cell.loc[:, "rna spots in pbody"] / (gene_Cell.loc[:, "nb_rna_out_nuc"] + gene_Cell.loc[:, "nb_rna_in_nuc"])).std()]

    #plot
    fig = gene_bar_plot(gene_list, mean_rna_per_pbody_list, title= title,path_output= path_output, ext=ext, show=show, close= close)




def RNApercentage_in_out_nucleus(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, plot_in_and_out_bars= True, path_output= None, show = True, close= True, ext= 'png', title = None) :

    join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
    if gene_list == None : gene_list = gdata.from_Acquisition_get_rna(Acquisition)

    std_list = []
    mean_value_inside = []
    mean_value_outside = []
    for gene in gene_list : 
        gene_Cell = join_frame[join_frame["rna name"] == gene]
        mean_value_inside += [gene_Cell.loc[:, "proportion_rna_in_nuc"].mean() * 100]
        mean_value_outside += [(1- gene_Cell.loc[:, "proportion_rna_in_nuc"].mean()) * 100]
        std_list += [gene_Cell.loc[:, "proportion_rna_in_nuc"].std() * 100]

    #plot
    ylabel = "Percentage of RNA found inside nucleus (%)"
    if plot_in_and_out_bars : fig = gene_bar_plot(gene_list, [mean_value_inside, mean_value_outside], [std_list]*2, legend = ["inside nuc", "outside nuc"])
    else : 
        fig = gene_bar_plot(gene_list, mean_value_inside, std_list, legend = "inside nuc", title= title, ylabel= ylabel, path_output= path_output, ext=ext, show=show, close= close)


################
# Violin plots #
################

def violin_rna_in_pbody(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, ext= 'png', title = None):
    """
    Work in progress #TODO
    """


    join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
    join_frame = join_frame.drop(axis= 0, index= join_frame[join_frame["pbody number"] == 0].index)
    
    if gene_list == None : 
        gene_list = gdata.from_Acquisition_get_rna(Acquisition)
    
    #Sort gene alphabetically
    gene_list = pd.DataFrame(columns = ["rna name"], data= gene_list).sort_values("rna name")
    gene_list = list(gene_list["rna name"])

    #Build array for plot
    violin_array = []
    for gene in gene_list:
        gene_frame = join_frame.query("`rna name` == '{0}'".format(gene))
        violin_array.append(gene_frame.loc[:, "rna spots in body"])
    print(violin_array)
    #plot
    fig = plt.figure(figsize= (100,10))
    ax = fig.gca()
    plt.xticks(range(len(gene_list)))
    xticks = ax.get_xticks()
    ax.set_xticks(xticks, labels= gene_list, rotation= 90)
    ax.axis(xmin=-0.5, xmax= len(gene_list)+0.5)
    fig.subplots_adjust(bottom= 2/fig.get_size_inches()[1])

    plt.violinplot(violin_array, positions= np.arange(len(gene_list)))
    if show : plt.show()




################################
## Malat1 / Dapi Signal plots ##
################################

def Malat_inNuc_asDapiIntensity(Cell: pd.DataFrame, projtype = 'MIP', summarize_type= 'median', out = False, plot_linear_regression= False,
                                path_output= None, show = True, ext= 'png', title = None) :
    """
    
    Scatter plot computed as Signal from DAPI channel (X-axis) VS Malat spots count in nucleus.
    Signal can be chosen from mean or median and computed from the Maximum Intensity Projection (MIP) or Mean projection.
    
    """



    #Select projection type from Cell Data
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype shoud either be 'mip' or 'mean', it is {0}.".format(projtype))

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    #Select malat spots in nucleus or in cytoplasm
    if out : Y = 'malat1 spots out nucleus'
    else : Y = 'malat1 spots in nucleus'

    #Fetching data
    Cell["SignalArea"] = Cell[X] * Cell["nuc_area"]
    data = Cell.loc[:, ["SignalArea", Y]].sort_values("SignalArea")

    #Plot
    fig = plt.figure(figsize=(20,10))
    X_values, Y_values = np.array(data.iloc[:,0]), np.array(data.iloc[:,1])
    plt.plot(X_values, Y_values, 'r.', label= "Experimental Data")

    if plot_linear_regression:
        slope, intercept = simple_linear_regression(X_values,Y_values)
        regression = slope* X_values + intercept
        plt.plot(X_values, regression, 'b', label= "Linear regression \n{0}x + {1}".format(slope,intercept))
        plt.legend()

    ax = fig.gca()
    ax.axis(xmin= 0, ymin= 0)
    if title != None : plt.title(title)
    plt.xlabel(X)
    plt.ylabel(Y)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()


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


###############
# Histogramms #
###############

def hist_dapi_signal(Cell, projtype= 'MIP', summarize_type = 'median', path_output= None, show = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    From the Maximum Intensity Projection (MIP) or Mean projection computes the histogram of integrated signal within cells (signal*nucleus area).
    Signal can be chosen from mean or median.

    Parameters
    ----------

        projtype : "MIP" or "MEAN"
        summarize_type : "median" or "mean"
        axis_boundaries : kwargs
            boundaries for x and y axes. Expected None or at least one the following ('xmin'=x, 'xmax'=X, ymin='y',ymax='Y')

    """
    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    Cell = update.from_nucleus_malat_proportion_compute_CellullarCycleGroup(Cell, 0.5)
    dapi_signal = Cell.loc[:, [X, "Cellular_cycle (malat proportion)"]]
    g1 = dapi_signal.query("`Cellular_cycle (malat proportion)` == 'g1'").index
    g2 = dapi_signal.query("`Cellular_cycle (malat proportion)` == 'g2'").index
    dapi_signal[X] = Cell.loc[:,X] * Cell.loc[:,"nuc_area"]
    histogram(dapi_signal.loc[g1,X], color= 'green', barlabel= 'g1', close= False, show= False)
    histogram(dapi_signal.loc[g2, X], color='red', barlabel= 'g2', reset= False, xlabel="Dapi signal (Nucleus area * {0})".format(X), ylabel= "count", path_output=path_output, show=show, ext=ext, title=title, bins=bins, **axis_boundaries)


def hist_dapi_density(Cell, projtype= 'MIP', summarize_type = 'median', path_output= None, show = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    From the Maximum Intensity Projection (MIP) or Mean projection computes the histogram of dapi density within cells (signal/nucleus area).
    Signal can be chosen from mean or median.

    Parameters
    ----------

        projtype : "MIP" or "MEAN"
        summarize_type : "median" or "mean"
        axis_boundaries : kwargs
            boundaries for x and y axes. Expected None or at least one the following ('xmin'=x, 'xmax'=X, ymin='y',ymax='Y')

    """
    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    dapi_signal = Cell.loc[:,X] / Cell.loc[:,"nuc_area"]
    histogram(dapi_signal, xlabel="Dapi signal ({0}/ Nucleus area)".format(X), ylabel= "count", path_output=path_output, show=show, ext=ext, title=title, bins=bins, **axis_boundaries)



def hist_malat_count(Cell, location= 'nucleus', path_output= None, show = True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    Histogram of malat spots detected in nucleus, cytoplasm or both.
    location : 'nucleus', 'cytoplasm' or 'cell' (cell = nuc + cytoplasm)
    projtype : "MIP" or "MEAN"
    """
    Cell = update.from_nucleus_malat_proportion_compute_CellullarCycleGroup(Cell, 0.5)
    if location.upper() == "NUCLEUS" : X = "malat1 spots in nucleus"
    elif location.upper() == "CYTOPLASM" or location.upper() == "CELL" : X = "malat1 spots in cytoplasm"
    else : raise ValueError("Incorrect value for location parameter. Should be one of the following : 'nucleus', 'cytoplasm' or 'cell'.")
    dapi_signal = Cell.loc[:, [X, "Cellular_cycle (malat proportion)"]]
    g1 = dapi_signal.query("`Cellular_cycle (malat proportion)` == 'g1'").index
    g2 = dapi_signal.query("`Cellular_cycle (malat proportion)` == 'g2'").index
    if location.upper() == "CELL" : dapi_signal += Cell.loc[:, "malat1 spots in nucleus"]
    
    histogram(dapi_signal.loc[:, X], color= 'green',xlabel=X, ylabel= "count", show=show, path_output=path_output, ext=ext, title=title, bins=bins, **axis_boundaries)


def hist_in_nuc_malat_proportion(Cell, path_output= None, show = True, close = True, ext= 'png', title: str = None, bins= 500, **axis_boundaries) :
    """
    Histogram of malat proportion detected inside nucleus.
    """
    proportion = Cell.loc[:, "malat1 spots in nucleus"] / (Cell.loc[:, "malat1 spots in nucleus"] + Cell.loc[:, "malat1 spots in cytoplasm"])
    histogram(proportion, xlabel="malat spots proportion in nucleus", ylabel= "count", path_output=path_output, show=show, close=close, ext=ext, title=title, bins=bins, **axis_boundaries)




################
# P-Body plots #
################

def bar_P_body_detect_inside_nucleus(Cell, path_output= None, show = True, close = True, ext= 'png', title: str = None) :
    """
    Plot a 2 bar graphs : one bar for number of pbodies detected inside nucleus and one for number of pbodies detected in cytoplasm.
    """

    Number_in_cytoplasm = Cell.loc[:, "count pbody in cytoplasm"].sum()
    std_in_cytoplasm = Cell.loc[:, "count pbody in cytoplasm"].std()
    Number_in_nucleus = Cell.loc[:, "count pbody in nucleus"].sum()
    std_in_nucleus = Cell.loc[:, "count pbody in nucleus"].std()

    ylabel = "Count"

    gene_bar_plot(["P-bodies in cytoplasm", "P-Bodies in nucleus"], [Number_in_cytoplasm, Number_in_nucleus], [std_in_cytoplasm, std_in_nucleus],ylabel= ylabel, path_output=path_output, show=show, close=close, ext=ext, title=title)