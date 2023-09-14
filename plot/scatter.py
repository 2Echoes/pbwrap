"""
This submodules groups all function related to scatter plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.update as update
import CustomPandasFramework.PBody_project.get as get
import matplotlib.pyplot as plt
from CustomPandasFramework.integrity.Errors import MissingColumnError
from ..quantification.CurveAnalysis import simple_linear_regression
from .utils import save_plot, get_colors_list, annotate_plot, get_markers_generator, hide_overlapping_annotations

def Malat_inNuc_asDapiIntensity(CellularCycle_view: pd.DataFrame, plot_linear_regression= False, path_output= None, show = True, close= True, reset= True, ext= 'png', title = None, xlabel=None, ylabel= 'malat1 spots in nucleus', **kargs):
    """
    Nuclei malat spot count VS Dapi Integrated Signal
    """
    cellcycle_view = CellularCycle_view.copy().dropna(subset=["IntegratedSignal", "count_in_nuc", "count_in_cyto"])

    rna_list = cellcycle_view.index.get_level_values(0).unique() 

    #axis labels
    if xlabel == None : xlabel = "Integrated Dapi Signal"
    if ylabel == None : ylabel = "Malat spot count"

    #colors
    if not "color" in kargs :
        auto_color = True 
        color_gen =  iter(get_colors_list(len(rna_list)))
    else : auto_color = False


    for rna in rna_list :
        if auto_color : kargs["color"] = next(color_gen)
        X = cellcycle_view.loc[rna, "IntegratedSignal"]
        Y = cellcycle_view.loc[rna, "count_in_nuc"]
        scatter(X=X, Y=Y, xlabel=xlabel, ylabel=ylabel, show=False, close=False, reset=reset, title=title, **kargs)
        plt.scatter(x=[],y=[], label= rna, **kargs)
        reset = False
    
    if plot_linear_regression :
        slope, intercept = simple_linear_regression(X = cellcycle_view.loc[:, "IntegratedSignal"], Y= cellcycle_view.loc[:, "count_in_nuc"] + cellcycle_view.loc[:, "count_in_cyto"])
        xmin = X.min()
        xmax = X.max()
        xrange = np.linspace(xmin,xmax, 100)
        plt.plot(xrange, slope*xrange + intercept, label= 'Linear regression : {0}x + {1}'.format(round(slope,10), round(intercept,2)))

    plt.legend()
    if type(path_output) != type(None) : save_plot(path_output=path_output,ext=ext)
    if show : plt.show()
    if close : plt.close()

def count_Malat_per_Cell(CellCellular_cycle, xlabel= "Mean", ylabel= "Standard deviation", title= "Malat spots detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    """
    df = CellCellular_cycle.copy().dropna(subset= ["count_in_nuc", "count_in_cyto"])
    df["total"] = df["count_in_nuc"] + df["count_in_cyto"]
    rna_list = CellCellular_cycle.index.get_level_values(0).unique() 

    #colors
    if not "color" in kargs :
        auto_color = True 
        color_gen =  iter(get_colors_list(len(rna_list)))
    else : auto_color = False

    for rna in rna_list :
        if auto_color : kargs["color"] = next(color_gen)
        
        df_rna = df.loc[rna,:].groupby(by= "AcquisitionId")["total"]
        Y = df_rna.std()
        X = df_rna.mean()

        scatter(X=X, Y=Y, xlabel=xlabel, ylabel=ylabel, show=False, close=False, reset=reset, title=title, **kargs)
        plt.scatter(x=[],y=[], label= rna, **kargs)
        reset = False


def count_rna_per_Cell(detection_view, xlabel= "Mean", ylabel= "Standard deviation", title= "Rna spots detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 point per fov
    """

    rna_list = detection_view.index.get_level_values(1).unique() 

    #axis labels
    if xlabel == None : xlabel = "Integrated Dapi Signal"
    if ylabel == None : ylabel = "Malat spot count"

    #colors
    if not "color" in kargs :
        auto_color = True 
        color_gen =  iter(get_colors_list(len(rna_list)))
    else : auto_color = False

    for rna in rna_list :
        if auto_color : kargs["color"] = next(color_gen)
        group = detection_view.loc[("rna", rna), ["AcquisitionId", "count"]].groupby("AcquisitionId")
        X = group["count"].mean()
        Y = group["count"].std()
        scatter(X=X, Y=Y, xlabel=xlabel, ylabel=ylabel, show=False, close=False, reset=reset, title=title, **kargs)
        plt.scatter(x=[],y=[], label= rna, **kargs)
        reset = False

    
    plt.legend()
    if type(path_output) != type(None) : save_plot(path_output=path_output,ext=ext)
    if show : plt.show()
    if close : plt.close()

def count_pbody_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= "Mean", ylabel= "Standard deviation", title= "P-bodies detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 point per fov. Obsolete.
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])

    Df_Acquisition = pd.merge(left= Join_Cell.loc[:,["rna name","AcquisitionId", "pbody number"]].groupby(["rna name","AcquisitionId"]).mean()["pbody number"], 
                              right=Join_Cell.loc[:,["rna name","AcquisitionId", "pbody number"]].groupby(["rna name","AcquisitionId"]).std()["pbody number"]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={'pbody number_x' : 'mean', 'pbody number_y' : 'std'})
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_frame = Join_Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
    color_list = pd.DataFrame(columns = ["color"], data = get_colors_list(gene_number))
    gene_frame = pd.concat([gene_frame, color_list], axis= 1).drop(0, axis= 1)
    Df_Acquisition = pd.merge(left= Df_Acquisition, right= gene_frame, how= 'left', left_on= "rna name", right_on= "rna name")

    scatter(X= Df_Acquisition["mean"], Y = Df_Acquisition["std"], xlabel= xlabel, ylabel= ylabel, color= Df_Acquisition["color"], label= list(Df_Acquisition["rna name"]), title=title, reset=reset, close=close, show=show, path_output=path_output, ext=ext, **kargs)




def dapi_signal(Cell: pd.DataFrame, projtype= 'mean', summarize_type= 'mean', integrated_signal = False,
                xlabel= "Mean", ylabel= "Standard deviation", title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    Integrated signal -> True : signal value is multiplied with nucleus area (nm)
    """
    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    
    if integrated_signal:
        Cell["Integrated signal ({0})".format(X)] = (Cell[X] * Cell["nucleus area (nm^2)"])
        X = "Integrated signal ({0})".format(X)

    if title == None : title = X

    Df_Acquisition = pd.merge(left= Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).mean()[X], 
                              right=Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).std()[X]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={ X+'_x' : 'mean', X+'_y' : 'std'})
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_frame = Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Cell.value_counts(subset="rna name"))
    color_list = pd.DataFrame(columns = ["color"], data = get_colors_list(gene_number))
    gene_frame = pd.concat([gene_frame, color_list], axis= 1).drop(0, axis= 1)
    Df_Acquisition = pd.merge(left= Df_Acquisition, right= gene_frame, how= 'left', left_on= "rna name", right_on= "rna name")

    scatter(X= Df_Acquisition["mean"], Y = Df_Acquisition["std"], xlabel= xlabel, ylabel= ylabel, color= Df_Acquisition["color"], label= list(Df_Acquisition["rna name"]), title=title, reset=reset, close=close, show=show, path_output=path_output, ext=ext, **kargs)

def DapiSignal_vs_CellNumber(Cell: pd.DataFrame, projtype= 'mean', summarize_type= 'mean', integrated_signal = False,
                             xlabel= "Cell Number", ylabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    #Projtype
    if projtype.upper() == 'MIP' : Y = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : Y = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : Y += "median_signal"
    elif summarize_type.upper() == 'MEAN' : Y += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    if integrated_signal:
        Cell["Integrated signal ({0})".format(Y)] = (Cell[Y] * Cell["nucleus area (nm^2)"])
        Y = "Integrated signal ({0})".format(Y)

    if ylabel == None : ylabel = Y

    mean_df = Cell.groupby("AcquisitionId")[Y].mean()
    count_df = Cell.value_counts(subset= "AcquisitionId").reset_index(drop=False).rename(columns={0 : 'cell number'})
    df = pd.merge(mean_df, count_df, on= "AcquisitionId").sort_values(Y)
    
    plot(df["cell number"], df[Y], xlabel= xlabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)


def G1G2_SpotsInPbody_Quantif(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type= 'rna',
                          xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")
    
    if reset : fig, _ = plt.subplots(nrows= 1, ncols= 3, figsize= (45, 15))
    else : fig = plt.gcf()
    if type(title) != type(None) : plt.suptitle(title, fontsize= 35, fontweight= 'bold')
    if type(xlabel) != type(None) : fig.text(0.5, 0.04, xlabel, ha='center', fontsize= 65)
    
    RNA_idx = Spots.query("spots_type == '{0}'".format(spots_type)).index
    Spots = Spots.copy().loc[RNA_idx,:]
    
    Cell_df = update.removeCellsWithoutSpotsInPodies(Cell, Spots)
    gene_list = list(Cell_df.groupby(['rna name'])["id"].count().sort_index().index)

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

    red_genes = get.get_GenesWithlessThanXCells(Cell_df.query("cellular_cycle == 'g1'"),100,'id')
    red_genes += get.get_GenesWithlessThanXCells(Cell_df.query("cellular_cycle == 'g2'"),100,'id')
    red_index = [gene_list.index(gene) for gene in red_genes]
    for idx in red_index :
        edgecolors_list[idx] = 'red'
        linewidths_list[idx] = 1.2


    kargs["color"] = color_list
    kargs["linewidths"] = linewidths_list
    kargs["edgecolors"] = edgecolors_list
    kargs["alpha"] = 0.7

    plt.subplot(1,3,1)
    title= "Mean {0} number in Pbodies".format(spots_type)
    G1G2_RnaNumberInPbody(Cell=Cell_df, Spots=Spots, legend=True, title=title, **kargs)
    plt.subplot(1,3,2)
    title = "Mean {0} proportion in Pbodies".format(spots_type)
    G1G2_RnaProportionInPbody(Cell=Cell_df, Spots=Spots, legend=True, title=title, **kargs)
    ax =plt.subplot(1,3,3)
    G1G2_CellNumber(Cell=Cell_df, legend=True, **kargs)
    handle1 = plt.scatter([],[], color= "white", linewidths= 1.5, edgecolor='black', label= 'Genes with more \nthan 100 cells computed')
    handle2 = plt.scatter([],[], color= "white", linewidths= 2, edgecolor='red', label= 'Genes with less \nthan 100 cells computed')
    handles = [handle1, handle2]
    labels = ['Genes with more than 100 cells computed', 'Genes with less than 100 cells computed']
    fig.legend(handles, labels, loc='upper left', prop={'size': 20})
    
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()

def G1G2_CytoSpotsInPbody_Quantif(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type= 'rna',
                          xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")
    
    if reset : fig, _ = plt.subplots(nrows= 1, ncols= 3, figsize= (45, 15))
    else : fig = plt.gcf()
    if type(title) != type(None) : plt.suptitle(title, fontsize= 35, fontweight= 'bold')
    if type(xlabel) != type(None) : fig.text(0.5, 0.04, xlabel, ha='center', fontsize= 65)
    
    RNA_idx = Spots.query("spots_type == '{0}'".format(spots_type)).index
    Spots = Spots.copy().loc[RNA_idx,:]
    
    Cell_df = update.removeCellsWithoutSpotsInPodies(Cell, Spots)
    gene_list = list(Cell_df.groupby(['rna name'])["id"].count().sort_index().index)

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

    red_genes = get.get_GenesWithlessThanXCells(Cell_df.query("cellular_cycle == 'g1'"),100,'id')
    red_genes += get.get_GenesWithlessThanXCells(Cell_df.query("cellular_cycle == 'g2'"),100,'id')
    red_index = [gene_list.index(gene) for gene in red_genes]
    for idx in red_index :
        edgecolors_list[idx] = 'red'
        linewidths_list[idx] = 1.2


    kargs["color"] = color_list
    kargs["linewidths"] = linewidths_list
    kargs["edgecolors"] = edgecolors_list
    kargs["alpha"] = 0.7

    plt.subplot(1,3,1)
    title= "Mean {0} number in Pbodies".format(spots_type)
    G1G2_cyto_spots_InPbody(Cell=Cell_df, Spots=Spots, legend=True, title=title, **kargs)
    plt.subplot(1,3,2)
    title = "Mean {0} proportion in Pbodies".format(spots_type)
    G1G2_cyto_spots_InPbody_proportion(Cell=Cell_df, Spots=Spots, legend=True, title=title, **kargs)
    ax =plt.subplot(1,3,3)
    G1G2_CellNumber(Cell=Cell_df, legend=True, **kargs)
    handle1 = plt.scatter([],[], color= "white", linewidths= 1.5, edgecolor='black', label= 'Genes with more \nthan 100 cells computed')
    handle2 = plt.scatter([],[], color= "white", linewidths= 2, edgecolor='red', label= 'Genes with less \nthan 100 cells computed')
    handles = [handle1, handle2]
    labels = ['Genes with more than 100 cells computed', 'Genes with less than 100 cells computed']
    fig.legend(handles, labels, loc='upper left', prop={'size': 20})
    
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()



def G1G2_Spots_Quantif(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type = 'rna',
                          xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")
    
    if reset : fig, _ = plt.subplots(nrows= 1, ncols= 2, figsize= (30, 15))
    else : fig = plt.gcf()
    if type(title) != type(None) : plt.suptitle(title, fontsize= 35, fontweight= 'bold')
    if type(xlabel) != type(None) : fig.text(0.5, 0.04, xlabel, ha='center', fontsize= 65)

    RNA_idx = Spots.query("spots_type == '{0}'".format(spots_type)).index
    Spots = Spots.copy().loc[RNA_idx,:]

    Cell_df = update.removeCellsWithoutSpotsInPodies(Cell, Spots)
    gene_list = list(Cell_df.groupby(['rna name'])["id"].count().sort_index().index)

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

    red_genes = get.get_GenesWithlessThanXCells(Cell_df.query("cellular_cycle == 'g1'"),100,'id')
    red_genes += get.get_GenesWithlessThanXCells(Cell_df.query("cellular_cycle == 'g2'"),100,'id')
    red_index = [gene_list.index(gene) for gene in red_genes]
    for idx in red_index :
        edgecolors_list[idx] = 'red'
        linewidths_list[idx] = 1.2


    kargs["color"] = color_list
    kargs["linewidths"] = linewidths_list
    kargs["edgecolors"] = edgecolors_list
    kargs["alpha"] = 0.7

    plt.subplot(1,2,1)
    title = "{0} spots per cell".format(spots_type)
    G1G2_rna_per_cell(Cell=Cell_df, Spots=Spots, legend=True, title=title, **kargs)
    plt.subplot(1,2,2)
    G1G2_CellNumber(Cell=Cell_df, legend=True, **kargs)
    handle1 = plt.scatter([],[], color= "white", linewidths= 1.5, edgecolor='black', label= 'Genes with more \nthan 100 cells computed')
    handle2 = plt.scatter([],[], color= "white", linewidths= 2, edgecolor='red', label= 'Genes with less \nthan 100 cells computed')
    handles = [handle1, handle2]
    labels = ['Genes with more than 100 cells computed', 'Genes with less than 100 cells computed']
    fig.legend(handles, labels, loc='upper left', prop={'size': 20})
    
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


def G1G2_RnaNumberInPbody(Cell: pd.DataFrame, Spots: pd.DataFrame, 
                          xlabel= "G1", ylabel= "G2", title= "Mean rna number in Pbodies",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :

    if 'rna name' not in Spots.columns : raise MissingColumnError("'rna name' column is missing from Spots DF : consider using update.AddRnaName")

    Spots_DF = pd.merge(Spots, Cell.loc[:,["id", "cellular_cycle"]], how= 'left', left_on= "CellId", right_on= 'id')
    Spots_DF = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["PbodyId"].count()
    drop_idx = Spots_DF[Spots_DF == 0].index
    Spots_DF = Spots_DF.drop(drop_idx)
    gene_list = Spots_DF.sort_index().index.get_level_values(0).unique()

    if reset : plt.figure(figsize=(20,20))

    kargs_copy = kargs.copy()
    del kargs_copy["color"],kargs_copy["linewidths"],kargs_copy["edgecolors"]


    markers_gen = get_markers_generator()
    annotation_list = []
    for gene, color, lw, edgecolor in zip(gene_list, kargs["color"], kargs["linewidths"], kargs["edgecolors"]) :
        marker = next(markers_gen)
        DF = Spots_DF.loc[gene,:]
        index_lvl0 = DF.index.get_level_values(0).unique()
        if "g1" in index_lvl0 : g1_mean = DF.loc["g1",:].mean()
        else : g1_mean = 0
        if "g2" in index_lvl0 : g2_mean = DF.loc["g2",:].mean()
        else : g2_mean = 0
        plt.scatter(x= g1_mean, y= g2_mean, color = color, label= gene, linewidths=lw, marker=marker, edgecolors= edgecolor, s= 60, **kargs_copy)
        annotation_list.append(plt.text(x= g1_mean*0.98, y = g2_mean*1.01, s= gene, size= 10))

    hide_overlapping_annotations(*annotation_list)
    if legend : plt.legend(ncols= 4)
    plt.axis("square")
    # plt.axis([0,300,0,300])
    if type(xlabel) != type(None) : plt.xlabel(xlabel)
    if type(ylabel) != type(None) : plt.ylabel(ylabel)
    if type(title) != type(None) : plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()



def G1G2_RnaProportionInPbody(Cell: pd.DataFrame, Spots: pd.DataFrame, 
                          xlabel= "G1", ylabel= "G2", title= "Mean rna proportion in Pbodies",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :

    if 'rna name' not in Spots.columns : raise MissingColumnError("'rna name' column is missing from Spots DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = pd.merge(Spots, Cell.loc[:,["id", "cellular_cycle"]], how= 'left', left_on= "CellId", right_on= 'id').drop('id_y', axis= 1).rename(columns= {'id_x' : 'id'})
    inPbody = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["PbodyId"].count()
    total = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["id"].count()
    # proportion = inPbody/total
    # drop_idx = proportion[proportion == 0].index
    # proportion = proportion.drop(drop_idx, axis=0)
    # assert all(proportion <= 1), "Error : proportion > 1 found while computing Spots proportion in Pbody"

    # Spots_DF = proportion
    # drop_idx = Spots_DF[Spots_DF == 0].index
    # Spots_DF = Spots_DF.drop(drop_idx)
    gene_list = total.sort_index().index.get_level_values(0).unique()

    if reset : plt.figure(figsize=(20,20))
    
    kargs_copy = kargs.copy()
    del kargs_copy["color"],kargs_copy["linewidths"],kargs_copy["edgecolors"]
    marker_gen = get_markers_generator()
    annotation_list = []
    for gene, color, lw, edgecolor in zip(gene_list, kargs["color"], kargs["linewidths"], kargs["edgecolors"]) :
        marker = next(marker_gen)
        DF_total = total.loc[gene,:]
        DF_inPbody= inPbody.loc[gene,:]
        DF =DF_total#Spots_DF.loc[gene,:]
        index_lvl0 = DF.index.get_level_values(0).unique()
        if "g1" in index_lvl0 : 
            g1_mean = DF_inPbody.loc["g1",:].sum() / DF_total.loc["g1",:].sum()
            if g1_mean == np.inf : print("Mean equals inf : \n",DF.loc["g1",:])
        else : g1_mean = 0
        if "g2" in index_lvl0 : 
            g2_mean = DF_inPbody.loc["g2",:].sum() / DF_total.loc["g2",:].sum()
        else : g2_mean = 0
        plt.scatter(x= g1_mean, y= g2_mean, color = color, label= gene, marker=marker, linewidths=lw, edgecolors= edgecolor, s=60, **kargs_copy)
        annotation_list.append(plt.text(x= g1_mean*0.98, y = g2_mean*1.02, s= gene, size= 7))

    hide_overlapping_annotations(*annotation_list)
    if legend : plt.legend(ncols= 4)
    axes = plt.axis('square')
    # plt.axis([0,0.6,0,0.6])
    if type(xlabel) != type(None) : plt.xlabel(xlabel)
    if type(ylabel) != type(None) : plt.ylabel(ylabel)
    if type(title) != type(None) : plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


def G1G2_CellNumber(Cell: pd.DataFrame,
                          xlabel= "G1", ylabel= "G2", title= "Number of cell computed", legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :

    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    cell_number = Cell.groupby(["rna name", "cellular_cycle"])["id"].count()
    gene_list = cell_number.sort_index().index.get_level_values(0).unique()

    if reset : plt.figure(figsize=(20,20))
 
    kargs_copy = kargs.copy()
    del kargs_copy["color"],kargs_copy["linewidths"],kargs_copy["edgecolors"]
    marker_gen = get_markers_generator()
    
    marker_gen = get_markers_generator()
    annotation_list = []
    for gene, color, lw, edgecolor in zip(gene_list, kargs["color"], kargs["linewidths"], kargs["edgecolors"]) :
        marker = next(marker_gen)
        DF = cell_number.loc[gene,:]
        index_lvl0 = DF.index
        if "g1" in index_lvl0 : g1_mean = DF.at["g1"]
        else : g1_mean = 0
        if "g2" in index_lvl0 : g2_mean = DF.at["g2"]
        else : g2_mean = 0

        plt.scatter(x= g1_mean, y= g2_mean, color = color, label= gene, marker=marker, linewidths=lw, edgecolors= edgecolor, **kargs_copy)
        annotation_list.append(plt.text(x= g1_mean*0.98, y = g2_mean*1.02, s= gene, size= 7))
    
    if legend : plt.legend(ncols= 4)
    hide_overlapping_annotations(*annotation_list)
    plt.axis('square')
    if type(xlabel) != type(None) : plt.xlabel(xlabel)
    if type(ylabel) != type(None) : plt.ylabel(ylabel)
    if type(title) != type(None) : plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


def G1G2_rna_per_cell(Cell: pd.DataFrame, Spots: pd.DataFrame, 
                          xlabel= "G1", ylabel= "G2", title= "Mean rna proportion in Pbodies",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = pd.merge(Spots, Cell.loc[:,["id", "cellular_cycle"]], how= 'left', left_on= "CellId", right_on= 'id').drop('id_y', axis= 1).rename(columns={'id_x' : 'id'})
    Spots_DF = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["id"].count()
    gene_list = Spots_DF.sort_index().index.get_level_values(0).unique()

    if reset : plt.figure(figsize=(20,20))

    kargs_copy = kargs.copy()
    del kargs_copy['color'], kargs_copy["edgecolors"], kargs_copy['linewidths']

    markers_gen = get_markers_generator()
    annotation_list = []
    for gene, color, lw, edgecolor in zip(gene_list, kargs["color"], kargs["linewidths"], kargs["edgecolors"]) :
        marker = next(markers_gen)
        DF = Spots_DF.loc[gene,:]
        index_lvl0 = DF.index.get_level_values(0).unique()
        if "g1" in index_lvl0 : g1_mean = DF.loc["g1",:].mean()
        else : g1_mean = 0
        if "g2" in index_lvl0 : g2_mean = DF.loc["g2",:].mean()
        else : g2_mean = 0
        plt.scatter(x= g1_mean, y= g2_mean, color = color, label= gene, linewidths=lw, marker=marker, edgecolors= edgecolor, s= 60, **kargs_copy)
        annotation_list.append(plt.text(x= g1_mean*0.98, y = g2_mean*1.01, s= gene, size= 10))
    
    if legend : plt.legend(ncols= 4)
    hide_overlapping_annotations(*annotation_list)
    plt.axis('square')
    # plt.axis([0,4000,0,4000])
    if type(xlabel) != type(None) : plt.xlabel(xlabel)
    if type(ylabel) != type(None) : plt.ylabel(ylabel)
    if type(title) != type(None) : plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()

def G1G2_cyto_spots_InPbody(Cell: pd.DataFrame, Spots: pd.DataFrame, 
                          xlabel= "G1", ylabel= "G2", title= "Mean rna proportion in Pbodies",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = pd.merge(Spots, Cell.loc[:,["id", "cellular_cycle"]], how= 'left', left_on= "CellId", right_on= 'id').drop('id_y', axis= 1).rename(columns={'id_x' : 'id'})
    Spots_DF["InCyto"] = 1 - ( Spots_DF["InNucleus"].astype(bool) | Spots_DF["PbodyId"].isna())
    Spots_DF = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["InCyto"].sum()
    gene_list = Spots_DF.sort_index().index.get_level_values(0).unique()

    if reset : plt.figure(figsize=(20,20))

    kargs_copy = kargs.copy()
    del kargs_copy['color'], kargs_copy["edgecolors"], kargs_copy['linewidths']

    markers_gen = get_markers_generator()
    annotation_list = []
    for gene, color, lw, edgecolor in zip(gene_list, kargs["color"], kargs["linewidths"], kargs["edgecolors"]) :
        marker = next(markers_gen)
        DF = Spots_DF.loc[gene,:]
        index_lvl0 = DF.index.get_level_values(0).unique()
        if "g1" in index_lvl0 : g1_mean = DF.loc["g1",:].mean()
        else : g1_mean = 0
        if "g2" in index_lvl0 : g2_mean = DF.loc["g2",:].mean()
        else : g2_mean = 0
        plt.scatter(x= g1_mean, y= g2_mean, color = color, label= gene, linewidths=lw, marker=marker, edgecolors= edgecolor, s= 60, **kargs_copy)
        annotation_list.append(plt.text(x= g1_mean*0.98, y = g2_mean*1.01, s= gene, size= 10))
    
    if legend : plt.legend(ncols= 4)
    hide_overlapping_annotations(*annotation_list)
    plt.axis('square')
    if type(xlabel) != type(None) : plt.xlabel(xlabel)
    if type(ylabel) != type(None) : plt.ylabel(ylabel)
    if type(title) != type(None) : plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


def G1G2_cyto_spots_InPbody_proportion(Cell: pd.DataFrame, Spots: pd.DataFrame, 
                          xlabel= "G1", ylabel= "G2", title= "Mean rna proportion in Pbodies",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = pd.merge(Spots, Cell.loc[:,["id", "cellular_cycle"]], how= 'left', left_on= "CellId", right_on= 'id').drop('id_y', axis= 1).rename(columns={'id_x' : 'id'})
    Spots_DF["InCyto"] = 1 - ( Spots_DF["InNucleus"].astype(bool) | Spots_DF["PbodyId"].isna())
    count = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["InCyto"].sum()
    total = Spots_DF.query("InNucleus == 0").groupby(["rna name", "cellular_cycle", "CellId"])["id"].count()
    Spots_DF = count/total
    assert Spots_DF[Spots_DF > 1].empty, "Error : proportion > 1 found while computing cytoplasmic spots proportion in Pbody"
    gene_list = Spots_DF.sort_index().index.get_level_values(0).unique()

    if reset : plt.figure(figsize=(20,20))

    kargs_copy = kargs.copy()
    del kargs_copy['color'], kargs_copy["edgecolors"], kargs_copy['linewidths']

    markers_gen = get_markers_generator()
    annotation_list = []
    for gene, color, lw, edgecolor in zip(gene_list, kargs["color"], kargs["linewidths"], kargs["edgecolors"]) :
        marker = next(markers_gen)
        DF = Spots_DF.loc[gene,:]
        index_lvl0 = DF.index.get_level_values(0).unique()
        if "g1" in index_lvl0 : g1_mean = DF.loc["g1",:].mean()
        else : g1_mean = 0
        if "g2" in index_lvl0 : g2_mean = DF.loc["g2",:].mean()
        else : g2_mean = 0
        plt.scatter(x= g1_mean, y= g2_mean, color = color, label= gene, linewidths=lw, marker=marker, edgecolors= edgecolor, s= 60, **kargs_copy)
        annotation_list.append(plt.text(x= g1_mean*0.98, y = g2_mean*1.01, s= gene, size= 10))
    
    if legend : plt.legend(ncols= 4)
    hide_overlapping_annotations(*annotation_list)
    plt.axis('square')
    if type(xlabel) != type(None) : plt.xlabel(xlabel)
    if type(ylabel) != type(None) : plt.ylabel(ylabel)
    if type(title) != type(None) : plt.title(title)

    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


## Base plot ##

def plot(X: np.ndarray, Y: np.ndarray, xlabel= None, ylabel= None, title= None, reset= False, close= False, show= True, path_output= None, ext ='png', **kargs) :
    #TODO does not function bc of label handling.
    """
    Default plot for points plots..query("`rna name` in ['NF1', 'PABPC1']")

    Parameters
    ----------
        data : sequence[float]
        
        **kargs :
            color

    """
    #auto kargs
    if "marker" not in kargs :
        kargs["marker"] = '.'
    if "ls" not in kargs :
        kargs['ls'] = ''


    if reset : fig = plt.figure(figsize=(20,10))
    else : fig = plt.gcf()

    plt.plot(X,Y, **kargs)
    if "label" in kargs :
        plt.legend()


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if path_output != None : save_plot(path_output=path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


    return fig



def scatter(X: np.ndarray, Y: np.ndarray, xlabel= None, ylabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    Default plot for scatter plots.

    Parameters
    ----------
        data : sequence[float]
        
        **kargs :
            color

    """

    if reset : fig = plt.figure(figsize=(20,20))
    else : fig = plt.gcf()


    #Auto kargs :
    if not 'edgecolors' in kargs :
        kargs['edgecolors'] = 'black'
    
    if "label" in kargs and "color" in kargs :
        label,color = kargs["label"], kargs["color"]
        del kargs["label"], kargs["color"]
        set_legend(labels= label, colors= color, **kargs)
        kargs["label"], kargs["color"] = label,color
        del label,color


    plt.scatter(X,Y, **kargs)
    
    plt.axis('tight')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    if path_output != None : save_plot(path_output=path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


    return fig


def set_legend(labels, colors, column_number = 3, loc = None, **kargs) :
    df = pd.DataFrame(data = {"label" : labels, "color" : colors})
    df = df.value_counts(subset=["label", "color"]).reset_index(drop= False)
    for label, color in zip(df["label"], df["color"]) :
        plt.scatter([],[],label = label, color = color, **kargs)
    

    plt.legend(ncol= column_number, loc=loc,)