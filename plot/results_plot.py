import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
from .utils import save_plot, gene_bar_plot, get_simple_linear_regression

def threshold(Acquisition: pd.DataFrame, rna_list:'list[str]' = None, path_output= None, show = True, ext= 'png', title = "threshold plot") :
    
    #Computing RNA mean threshold and var :
    if rna_list == None : rna_list =  gdata.from_Acquisition_get_rna(Acquisition)
    elif type(rna_list) == str : rna_list = [rna_list]

    threshold_list = []
    std_list = []
    for rna in rna_list :
        threshold_list += [Acquisition[Acquisition["rna name"] == rna].loc[:,"RNA spot threshold"].mean()]
        std_list += [Acquisition[Acquisition["rna name"] == rna].loc[:,"RNA spot threshold"].std()]

    
    fig = gene_bar_plot(rna_list, threshold_list, errors= std_list)
    plt.title(title)
    """
    #Plot
    color_list = ['red','blue','green','orange','purple','brown','cyan'] * (round(len(rna_list)/7) + 1)
    fig = plt.figure(figsize= (20,10))
    plt.bar(rna_list, threshold_list, yerr= std_list, capsize= 3, color= color_list[:len(rna_list)], width= 1)
    plt.axis(ymin= 0)
    
    
    #
    plt.xticks(range(len(rna_list)))
    ax = fig.gca()
    xticks = ax.get_xticks()
    ax.set_xticks(xticks, labels= rna_list, rotation= 90)
    ax.axis(xmin=-0.5, xmax= len(rna_list)+0.5, ymin= 0)
    fig.subplots_adjust(bottom= 2/fig.get_size_inches()[1])"""
    
    if path_output != None :
        save_plot(path_output= path_output, ext= ext)

    if show : plt.show()
    plt.close()



def spots_per_cell(Acquisition: pd.DataFrame, Cell: pd.DataFrame, spot_type = 'rna', rna_list: 'list[str]' = None, path_output= None, show = True, ext= 'png', title = "RNA per cell") :
    
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
    fig = gene_bar_plot(rna_list, threshold_list, std_list)
    plt.title(title)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()


def cluster_per_cell() :
    pass #no cluster data in frame

def RNA_in_pbody(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, ext= 'png', title = "RNA in pbody"):

    join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
    join_frame = join_frame.drop(axis= 0, index= join_frame[join_frame["pbody number"] == 0].index)
    if gene_list == None : gene_list = gdata.from_Acquisition_get_rna(Acquisition)
    
    std_list= []
    mean_rna_per_pbody_list = []
    for gene in gene_list : 
        gene_Cell = join_frame [join_frame["rna name"] == gene]
        mean_rna_per_pbody_list += [(gene_Cell.loc[:, "rna spots in body"] / gene_Cell.loc[:, "pbody number"]).mean()]
        std_list += [(gene_Cell.loc[:, "rna spots in body"] / gene_Cell.loc[:, "pbody number"]).std()]

    #plot
    fig = gene_bar_plot(gene_list, mean_rna_per_pbody_list, std_list)
    plt.title(title)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()

def RNApercentage_in_out_nucleus(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, ext= 'png', title = "Proportion of RNA inside and outside nucleus") :

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
    fig = gene_bar_plot(gene_list, [mean_value_inside, mean_value_outside], [std_list]*2, legend = ["inside nuc", "outside nuc"])
    plt.title(title)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()


def Malat_inNuc_asDapiIntensity(Cell: pd.DataFrame, projtype = 'MIP', out = False, plot_linear_regression= False,
                                path_output= None, show = True, ext= 'png', title = None) :

    #Select right projection type from Cell Data
    if projtype.upper() == 'MIP' : X = "Mean Intensity (MIP)"
    elif projtype.upper() == 'MEAN' : X = "Mean Intensity (MeanProj)"
    else : raise ValueError("projtype shoud either be 'mip' or 'mean', it is {0}.".format(projtype))

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
        slope, intercept = get_simple_linear_regression(X_values,Y_values)
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


def hist_dapi_signal(Cell, projtype= 'MIP', path_output= None, show = True, ext= 'png', title: str = None) :
    """
    projtype : "MIP" or "MEAN"
    """
    if projtype.upper() == 'MIP' : X = "Mean Intensity (MIP)"
    elif projtype.upper() == 'MEAN' : X = "Mean Intensity (MeanProj)"
    dapi_signal = Cell.loc[:,X] * Cell.loc[:,"nuc_area"]

    fig = plt.figure(figsize= (20,10))
    hist = plt.hist(dapi_signal, bins= 1000)
    plt.xlabel("Dapi signal (Nucleus area * mean signal)")
    plt.ylabel("Count")
    if title != None : plt.title(title)
    
    
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()



def hist_malat_count(Cell, out_nucleus= False, path_output= None, show = True, ext= 'png', title: str = None) :
    """
    projtype : "MIP" or "MEAN"
    """
    if out_nucleus : X = "malat1 spots in cytoplasm"
    else : X = "malat1 spots in nucleus"
    dapi_signal = Cell.loc[:, X]
    fig = plt.figure(figsize= (20,10))
    hist = plt.hist(dapi_signal, bins= 1000)
    plt.xlabel(X)
    plt.ylabel("Count")
    if title != None : plt.title(title)
    
    
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()



def hist_in_nuc_malat_proportion(Cell, path_output= None, show = True, ext= 'png', title: str = None) :
    """
    projtype : "MIP" or "MEAN"
    """
    proportion = Cell.loc[:, "malat1 spots in nucleus"] / (Cell.loc[:, "malat1 spots in nucleus"] + Cell.loc[:, "malat1 spots in cytoplasm"])
    fig = plt.figure(figsize= (20,10))
    hist = plt.hist(proportion * 100, bins= 1000)
    plt.xlabel("malat spots proportion in nucleus")
    plt.ylabel("Count")
    if title != None : plt.title(title)
    
    
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()