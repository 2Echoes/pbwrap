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

    #Plot
    color_list = ['red','blue','green','orange','purple','brown','cyan'] * (round(len(rna_list)/7) + 1)
    fig = plt.figure(figsize= (20,10))
    plt.bar(rna_list, threshold_list, yerr= std_list, capsize= 3, color= color_list[:len(rna_list)], width= 1)
    plt.axis(ymin= 0)
    plt.title(title)
    
    #
    plt.xticks(range(len(rna_list)))
    ax = fig.gca()
    xticks = ax.get_xticks()
    ax.set_xticks(xticks, labels= rna_list, rotation= 90)
    ax.axis(xmin=-0.5, xmax= len(rna_list)+0.5, ymin= 0)
    fig.subplots_adjust(bottom= 2/fig.get_size_inches()[1])
    
    if path_output != None :
        save_plot(path_output= path_output, ext= ext)

    if show : plt.show()
    plt.close()



def rna_per_cell(Acquisition: pd.DataFrame, Cell: pd.DataFrame, rna_list: 'list[str]' = None, path_output= None, show = True, ext= 'png', title = "RNA per cell") :
    
    
    if rna_list == None : rna_list =  gdata.from_Acquisition_get_rna(Acquisition)
    elif type(rna_list) == str : rna_list = [rna_list]

    Cell = gdata.from_rna_get_Cells(rna= rna_list, Cell= Cell, Acquisition= Acquisition)
    Cell["rna number"] = Cell["nb_rna_out_nuc"] + Cell["nb_rna_in_nuc"]

    #Computing mean values and std
    threshold_list = []
    std_list = []
    for rna in rna_list :
        threshold_list += [Cell[Cell["rna name"] == rna].loc[:,"rna number"].mean()]
        std_list += [Cell[Cell["rna name"] == rna].loc[:,"rna number"].std()]

    #plot
    fig = gene_bar_plot(rna_list, threshold_list, std_list)
    plt.title(title)
    if path_output != None : save_plot(path_output, ext)
    if show : plt.show()
    plt.close()


def RNA_in_pody():
    pass