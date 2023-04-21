import pandas as pd
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
from .utils import save_plot, gene_bar_plot

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

def RNApercentage_in_out_nlucleus(Acquisition: pd.DataFrame, Cell: pd.DataFrame, gene_list: 'list[str]' = None, path_output= None, show = True, ext= 'png', title = "Proportion on RNA inside or outside nucleus") :

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
    fig = gene_bar_plot(gene_list, [mean_value_inside, mean_value_outside], [std_list,std_list])
    plt.title(title)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()