"""
This submodules groups all function related to violin plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import pbwrap.data.getdata as gdata
import matplotlib.pyplot as plt


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
