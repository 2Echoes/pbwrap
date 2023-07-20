"""
This submodules groups all function related to scatter plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import CustomPandasFramework.PBody_project.update as update
import matplotlib.pyplot as plt
from ..quantification.CurveAnalysis import simple_linear_regression
from .utils import save_plot, get_colors_list, identity
from .decorators import plot_curve

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


def count_rna_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= "Mean", ylabel= "Standard deviation", title= "Rna spots detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
    """
    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    Join_Cell["count_rna_per_cell"] = (Join_Cell["nb_rna_out_nuc"] + Join_Cell["nb_rna_in_nuc"])

    Df_Acquisition = pd.merge(left= Join_Cell.loc[:,["rna name","AcquisitionId", "count_rna_per_cell"]].groupby(["rna name","AcquisitionId"]).mean()["count_rna_per_cell"], 
                              right=Join_Cell.loc[:,["rna name","AcquisitionId", "count_rna_per_cell"]].groupby(["rna name","AcquisitionId"]).std()["count_rna_per_cell"]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={'count_rna_per_cell_x' : 'mean', 'count_rna_per_cell_y' : 'std'})
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_frame = Join_Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
    color_list = pd.DataFrame(columns = ["color"], data = get_colors_list(gene_number))
    gene_frame = pd.concat([gene_frame, color_list], axis= 1).drop(0, axis= 1)
    Df_Acquisition = pd.merge(left= Df_Acquisition, right= gene_frame, how= 'left', left_on= "rna name", right_on= "rna name")

    scatter(X= Df_Acquisition["mean"], Y = Df_Acquisition["std"], xlabel= xlabel, ylabel= ylabel, color= Df_Acquisition["color"], label= list(Df_Acquisition["rna name"]), title=title, reset=reset, close=close, show=show, path_output=path_output, ext=ext, **kargs)


def count_pbody_per_Cell(Cell: pd.DataFrame, Acquisition: pd.DataFrame, xlabel= "Mean", ylabel= "Standard deviation", title= "P-bodies detected per Fov", reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    1 box per gene
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




def dapi_signal(Cell: pd.DataFrame, Acquisition: pd.DataFrame, projtype= 'mean', summarize_type= 'mean', integrated_signal = False,
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

    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    
    if integrated_signal:
        Join_Cell["Integrated signal ({0})".format(X)] = (Join_Cell[X] * Join_Cell["nucleus area (nm^2)"])
        X = "Integrated signal ({0})".format(X)

    if title == None : title = X

    Df_Acquisition = pd.merge(left= Join_Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).mean()[X], 
                              right=Join_Cell.loc[:,["rna name","AcquisitionId", X]].groupby(["rna name","AcquisitionId"]).std()[X]
                              ,left_on= ('rna name',"AcquisitionId"), right_on= ('rna name', "AcquisitionId")
                              ).rename(columns={ X+'_x' : 'mean', X+'_y' : 'std'})
    Df_Acquisition = Df_Acquisition.reset_index(drop= False).sort_values("rna name")
    gene_frame = Join_Cell.value_counts(subset="rna name").reset_index(drop= False)
    gene_number = len(Join_Cell.value_counts(subset="rna name"))
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
    
    plot(df["cell number"], df[Y], xlabel= xlabel, ylabel= ylabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)


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












    ################### OBSOLETTE ######################
"""



def Malat_inNuc_asDapiIntensity(Acquisition: pd.DataFrame, Cell: pd.DataFrame, Spots, projtype = 'MIP', summarize_type= 'median', integrated_signal=True, plot_linear_regression= False,
                                path_output= None, show = True, close= True, reset= True, ext= 'png', title = None, xlabel=None, ylabel= 'malat1 spots in nucleus', **kargs) :

    #Projtype
    if projtype.upper() == 'MIP' : X = "nucleus_mip_"
    elif projtype.upper() == 'MEAN' : X = "nucleus_mean_"
    else : raise ValueError("projtype should either be 'mip' or 'mean'.")

    #Summarize type
    if summarize_type.upper() == 'MEDIAN' : X += "median_signal"
    elif summarize_type.upper() == 'MEAN' : X += "mean_signal"
    else : raise ValueError("summarize_type should either be 'median' or 'mean'.")

    Join_Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    
    if integrated_signal:
        Join_Cell["integrated signal"] = (Join_Cell[X] * Join_Cell["nucleus area (nm^2)"])
        X = "integrated signal"
    
    if xlabel == None : xlabel = X

    rna_list = Join_Cell.value_counts("rna name")

    if not "color" in kargs : 
        color_df = pd.DataFrame({'rna name' : rna_list.index , 'color' : get_colors_list(len(rna_list))})
        Join_Cell = pd.merge(left= Join_Cell, right= color_df, how= 'left', on= 'rna name')
        kargs["color"] = Join_Cell["color"]

    
    kargs["label"] = Join_Cell["rna name"]

    if plot_linear_regression :
        
        if reset : 
            plt.figure(figsize=(20,10))
            reset = False

        slope, intercept = simple_linear_regression(X = Join_Cell[X], Y= Join_Cell['malat1 spots in nucleus'])
        xmin = Join_Cell[X].min()
        xmax = Join_Cell[X].max()
        xrange = np.linspace(xmin,xmax, 100)
        plt.plot(xrange, slope*xrange + intercept, label= 'Linear regression : {0}x + {1}'.format(round(slope,10), round(intercept,2)))


    scatter(X= Join_Cell[X], Y= Join_Cell['malat1 spots in nucleus'], xlabel=xlabel, ylabel=ylabel, show=show, close=close, reset=reset, path_output=path_output, ext=ext, title=title, **kargs)
"""