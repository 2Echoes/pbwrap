"""
This submodules groups all function related to scatter plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import warnings
import CustomPandasFramework.PBody_project.update as update
import CustomPandasFramework.PBody_project.get as get
import matplotlib.pyplot as plt
from CustomPandasFramework.integrity.Errors import MissingColumnError
from ..quantification.CurveAnalysis import simple_linear_regression
from .utils import save_plot, get_colors_list, annotate_plot, get_markers_generator, hide_overlapping_annotations
from .g1g2_layouts import _Layout_Quantif_plots, _G1G2_main_legend_layout, G1G2_plot


def RnaCount_VS_CellSize(Cell: pd.DataFrame, Spots: pd.DataFrame, plot_linear_regression= False, path_output= None, show = True, close= True, reset= True, title = None) :
    
    if 'rna name' not in Cell.columns : raise KeyError("rna name wasn't find in Cell columns")
    if 'rna name' not in Spots.columns : raise KeyError("rna name wasn't find in Spots columns")

    #data sorting
    malat1_idx = Spots.query("spots_type == 'rna'").index
    Spots_df = Spots.loc[malat1_idx, ['rna name', 'id', 'CellId']]

    signal_df = Cell.loc[:, ["rna name", "id", "cell_area"]].rename(columns={'id' : 'CellId'}).set_index(['rna name', 'CellId']).sort_index()

    rna_grouper = Spots_df.groupby(["rna name", "CellId"])
    rna_total = rna_grouper['id'].count()
    rna_total = rna_total.rename('rna_count')
    gene_list = rna_total.index.get_level_values(0).unique()

    #plot layout
    plot_number = len(gene_list)
    root = np.sqrt(plot_number)
    if root - np.floor(root) == 0 :
        n_lin = n_col = root
    elif root - np.floor(root) < 0.5 :
        n_lin = np.floor(root)
        n_col = np.floor(root) + 1
    else :
        n_lin = np.floor(root) + 1
        n_col = np.floor(root) + 1
    plt.close()
    if reset : fig, ax = plt.subplots(nrows= int(n_lin), ncols= int(n_col), figsize= (12*n_col, 12*n_lin))
    plot_idx=1
    color_list = iter(get_colors_list(plot_number))

    for gene in gene_list:
        color = next(color_list)
        plt.subplot(n_lin, n_col, plot_idx)
        x = signal_df[signal_df.index.isin(rna_total.index)].loc[gene,:]
        y = rna_total.loc[gene,:]
        plt.scatter(x = x, y= y, c= color, label= 'experimental data')
        plt.title(str(gene))
        plt.xlabel('Cell area')
        plt.ylabel('rna count')
        if plot_linear_regression :
            slope, intercept = simple_linear_regression(X = x, Y= y)
            xmin = signal_df.loc[gene,:].min()
            xmax = signal_df.loc[gene,:].max()
            xrange = np.linspace(xmin,xmax, 100)
            plt.plot(xrange, slope*xrange + intercept, label= 'Linear regression : {0}x + {1}'.format(round(slope,10), round(intercept,2)))

        plt.legend()
        xmin,xmax,ymin,ymax = plt.axis('tight')
        plt.axis([0,20000,0,400])
        plot_idx +=1

    if type(title) != type(None) : plt.suptitle(title, fontsize= 80, fontweight= 'bold')
    plt.tight_layout()
    if type(path_output) != type(None) : plt.savefig(path_output)
    if show : plt.show()
    if close : plt.close()




def Malat_inNuc_asDapiIntensity(Cell: pd.DataFrame, Spots: pd.DataFrame, plot_linear_regression= False, path_output= None, show = True, close= True, reset= True, title = None):
    """
    Nuclei malat spot count VS Dapi Integrated Signal.
    Cell without spots are filtered from plots.
    """
    if 'rna name' not in Cell.columns : raise KeyError("rna name wasn't find in Cell columns")
    if 'rna name' not in Spots.columns : raise KeyError("rna name wasn't find in Spots columns")

    #data sorting
    malat1_idx = Spots.query("spots_type == 'malat1'").index
    Spots_df = Spots.loc[malat1_idx, ['rna name', 'id', 'CellId', 'InNucleus']]
    if "IntegratedSignal" not in Cell.columns : signal_df = update.compute_IntegratedSignal(Cell).loc[:, ['rna name', "id", "IntegratedSignal"]].rename(columns={'id' : 'CellId'}).set_index(['rna name', 'CellId']).sort_index()
    else : signal_df = Cell.loc[:, ["rna name", "id", "IntegratedSignal"]].rename(columns={'id' : 'CellId'}).set_index(['rna name', 'CellId']).sort_index()

    malat1_grouper = Spots_df.groupby(["rna name", "CellId"])
    malat1_total = malat1_grouper['id'].count()
    malat1_count = (malat1_total).rename('InNucleus_malat1_count')
    gene_list = malat1_count.index.get_level_values(0).unique()

    #plot layout
    plot_number = len(gene_list)
    root = np.sqrt(plot_number)
    if root - np.floor(root) == 0 :
        n_lin = n_col = root
    elif root - np.floor(root) < 0.5 :
        n_lin = np.floor(root)
        n_col = np.floor(root) + 1
    else :
        n_lin = np.floor(root) + 1
        n_col = np.floor(root) + 1
    plt.close()
    if reset : fig, ax = plt.subplots(nrows= int(n_lin), ncols= int(n_col), figsize= (12*n_col, 12*n_lin))
    plot_idx=1
    color_list = iter(get_colors_list(plot_number))

    for gene in gene_list:
        color = next(color_list)
        plt.subplot(n_lin, n_col, plot_idx)
        x = signal_df[signal_df.index.isin(malat1_count.index)].loc[gene,:]
        y = malat1_count.loc[gene,:]
        plt.scatter(x = x, y= y, c= color, label= 'experimental data')
        plt.title(str(gene))
        plt.xlabel('Mean integrated signal')
        plt.ylabel('Malat1 count')
        if plot_linear_regression :
            slope, intercept = simple_linear_regression(X = x, Y= y)
            xmin = signal_df.loc[gene,:].min()
            xmax = signal_df.loc[gene,:].max()
            xrange = np.linspace(xmin,xmax, 100)
            plt.plot(xrange, slope*xrange + intercept, label= 'Linear regression : {0}x + {1}'.format(round(slope,10), round(intercept,2)))

        plt.legend()
        xmin,xmax,ymin,ymax = plt.axis('tight')
        plt.axis([0,xmax,0,ymax])
        plot_idx +=1

    if type(title) != type(None) : plt.suptitle(title, fontsize= 80, fontweight= 'bold')
    plt.tight_layout()
    if type(path_output) != type(None) : plt.savefig(path_output)
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


def G1G2_SpotsInPbody_Quantif(Cell: pd.DataFrame, Pbody: pd.DataFrame, Spots: pd.DataFrame, spots_type= 'rna', distance_SpotPbody= 0,  cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    """
    Plot figure containing 3 G1/G2 plots : RnaNumberInPbody, RnaProportionInPbody, ComputuedCellNumber

    Parameters
    ----------
        Cell : pd.DataFrame
            Cell results dataframe.
        Pbody : pd.DataFrame
            Pbody results dataFrame.
        Spots : pd.DataFrame
            Spots results dataFrame
        spots_type : str
            type of spots considered during ploting either 'rna' or 'malat1'.
        distance_SpotPbody : int,float-like
            Distance up to which a spot can be considered in the closest P-body. 
            0 can be passed to get spots only within the mask.
            To know which distance settings have been computed during analysis use `get.get_pbody_spot_distance_parameters()`
        
    """

    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")
    if distance_SpotPbody not in get.get_pbody_spot_distance_parameters(Pbody) : raise ValueError("value passed for 'distance_SpotPbody' was not computed during analysis. Please choose one amongst {0}".format(get.get_pbody_spot_distance_parameters(Pbody)))
    
    #Data management init
    # SpotsType_idx = Spots.query("spots_type == '{0}'".format(spots_type)).index
    # Spots_df = Spots.loc[SpotsType_idx,:]
    Spots_df = Spots[Spots['spots_type'] == spots_type]
    
    spot_measure = '{0} {1}nm count'.format(spots_type, distance_SpotPbody)
    Cell_df = Cell[Cell.id.isin(Pbody.CellId)]
    gene_list = list(Cell_df.value_counts("rna name").index)
    gene_outlier_dict = {
        'Df' : Cell_df,
        'number' : 100,
        'pk' : 'id'
        }
    fig, _, kargs = _Layout_Quantif_plots(gene_list, gene_outlier_dict,
                                            xlabel= xlabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)

    #Subplots
    plt.subplot(1,3,1)
    title= "Mean {0} number in Pbodies".format(spots_type)
    G1G2_RnaNumberInPbody(Cell=Cell_df, Pbody=Pbody, spot_measure=spot_measure,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, title=title, **kargs)
    plt.subplot(1,3,2)
    title = "Mean {0} proportion in Pbodies".format(spots_type)
    G1G2_RnaProportionInPbody(Cell=Cell_df, Pbody=Pbody, Spots=Spots_df, spot_measure=spot_measure,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, title=title, **kargs)
    plt.subplot(1,3,3)
    G1G2_PbodiesPerCell(Cell=Cell_df, Pbody=Pbody,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, **kargs)

    #Main plot legend
    if type(gene_outlier_dict) != type(None) : fig = _G1G2_main_legend_layout(fig)
    
    #Output
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()

def G1G2_CytoSpotsInPbody_Quantif(Cell: pd.DataFrame, Spots: pd.DataFrame, Pbody:pd.DataFrame, spots_type= 'rna', cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    """
    Plot figure containing 3 G1/G2 plots : G1G2_cyto_spots_InPbody, G1G2_cyto_spots_InPbody_proportion, G1G2_CellNumber
    For this plot only spots strictly inside the pbodies labels are considered.

    Parameters
    ----------
        Cell : pd.DataFrame
            Cell results dataframe.
        Pbody : pd.DataFrameFalse
        spots_type : str
            type of spots plotted : 'rna' or 'malat1'
        
    """
    
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")
    
    RNA_idx = Spots.query("spots_type == '{0}'".format(spots_type)).index
    Spots_df = Spots.loc[RNA_idx,:]
    
    Cell_df = update.removeCellsWithoutSpotsInPodies(Cell, Spots) #Ok make sense
    Cell_df = Cell
    gene_list = list(Cell_df.groupby(['rna name'])["id"].count().sort_index().index)
    gene_outlier_dict = {
        'Df' : Cell_df,
        'number' : 100,
        'pk' : 'id'
        }

    fig,_,kargs = _Layout_Quantif_plots(gene_list, gene_outlier_dict,
                             xlabel= xlabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)

    plt.subplot(1,3,1)
    title= "Mean {0} number in Pbodies".format(spots_type)
    G1G2_cyto_spots_InPbody(Cell=Cell_df, Spots=Spots_df, spots_type=spots_type,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, title=title, **kargs)
    plt.subplot(1,3,2)
    title = "Mean {0} proportion in Pbodies".format(spots_type)
    G1G2_cyto_spots_InPbody_proportion(Cell=Cell_df, Spots=Spots_df, spots_type=spots_type,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, title=title, **kargs)
    plt.subplot(1,3,3)
    G1G2_PbodiesPerCell(Cell=Cell_df, Pbody=Pbody,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, **kargs)

    fig = _G1G2_main_legend_layout(fig)
    
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()



def G1G2_Spots_Quantif(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type = 'rna', cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    """
    Plot figure containing 3 G1/G2 plots : G1G2_spots_per_cell, total_spots_per_gene, G1G2_CellNumber

    Parameters
    ----------
        Cell : pd.DataFrame
            Cell results dataframe.
        Pbody : pd.DataFrameFalse
        spots_type : str
            type of spots plotted : 'rna' or 'malat1'
        
    """
    
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    gene_list = list(Cell.groupby(['rna name'])["id"].count().sort_index().index)
    gene_outlier_dict = {
        'Df' : Cell,
        'number' : 100,
        'pk' : 'id'
        }
    fig,_,kargs = _Layout_Quantif_plots(gene_list, gene_outlier_dict,
                             xlabel= xlabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)

    plt.subplot(1,3,1)
    title = "{0} spots per cell".format(spots_type)
    G1G2_spots_per_cell(Cell=Cell, Spots=Spots, spots_type= spots_type,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, **kargs)
    plt.subplot(1,3,2)
    G1G2_total_spotnumber(Cell, Spots, spots_type=spots_type,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend= True, **kargs)
    plt.subplot(1,3,3)
    G1G2_CellNumber(Cell=Cell, legend=True,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, **kargs)
    fig = _G1G2_main_legend_layout(fig)
    
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()

    
def G1G2_KIF1C_plateQuantif(Spots_dataframe: pd.DataFrame, spots_type : str, cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                            xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    """
    Here plate name is equivalent to rna name (code speaking) in other plots.
    """
    if 'plate name' not in Spots_dataframe : raise KeyError("Spots dataframe is missing `plate name`column.")
    if 'rna name' in Spots_dataframe : 
        warnings.warn("rna name column was found within DataFrame, this is not supposed to happened when computing KIF1C variability. It has been dropped but you shoud investigate if behavior is as expected.")
        Spots_dataframe = Spots_dataframe.drop("rna name", axis= 1)
    type_idx = Spots_dataframe.query('spots_type == "{0}"'.format(spots_type)).index
    Spots_frame = Spots_dataframe.loc[type_idx,:]
    Spots_frame = Spots_frame.rename(columns={"plate name" : 'rna name'})
    Spots_Series: pd.Series = Spots_frame.groupby(['rna name', 'cellular_cycle', 'CellId'])["id"].count()
    gene_outlier_dict = {
        'Df' : Spots_frame,
        'number' : 100,
        'pk' : 'CellId'
        }
    fig,_,kargs = _Layout_Quantif_plots(list(Spots_Series.index.get_level_values(0).unique()), gene_outlier_dict,
                             xlabel= xlabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)

    plt.subplot(1,3,1)
    title = "KIF1C spots total"
    serie = Spots_Series.groupby(level=[0,1]).count()
    G1G2_plot(serie, cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, title=title)
    plt.subplot(1,3,2)
    title = "KIF1C mean spots number per cell"
    serie = Spots_Series
    G1G2_plot(serie, cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,title=title)
    plt.subplot(1,3,3)
    title = "Number of cells computed"
    serie = Spots_Series.groupby(level=[0,1]).count()
    G1G2_plot(serie, cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, title=title)

    fig = _G1G2_main_legend_layout(fig)
    
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()


def G1G2_spot_area_Quantif(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type : str, cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                            xlabel= None, title= None, reset= True, close= True, show= True, path_output= None, ext ='png', **kargs) :
    
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    gene_list = list(Cell.groupby(['rna name'])["id"].count().sort_index().index)
    gene_outlier_dict = {
        'Df' : Cell,
        'number' : 100,
        'pk' : 'id'
        }
    fig,_,kargs = _Layout_Quantif_plots(gene_list, gene_outlier_dict,
                             xlabel= xlabel, title= title, reset= reset, close= close, show= show, path_output= path_output, ext =ext, **kargs)

    plt.subplot(1,3,1)
    title = "{0} spots per cell".format(spots_type)
    G1G2_spots_per_cell(Cell=Cell, Spots=Spots, spots_type= spots_type,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend=True, **kargs)
    plt.subplot(1,3,2)
    G1G2_spots_per_area(Cell, Spots, spots_type=spots_type,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, legend= True, **kargs)
    plt.subplot(1,3,3)
    G1G2_CellNumber(Cell=Cell, legend=True,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y, **kargs)
    fig = _G1G2_main_legend_layout(fig)
    
    if type(path_output) != type(None) : save_plot(path_output, ext=ext)
    if show : plt.show()
    if close : plt.close()
    pass


def G1G2_RnaNumberInPbody(Cell: pd.DataFrame, Pbody: pd.DataFrame, spot_measure, cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          title= "Mean rna number in Pbodies per cell",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :

    if 'rna name' not in Pbody.columns : raise MissingColumnError("'rna name' column is missing from Spots DF : consider using update.AddRnaName")

    Pbody_DF = update.AddCellularCycle(Pbody, Cell)
    Pbody_DF = Pbody_DF.groupby(["rna name", "cellular_cycle", "CellId"])[spot_measure].sum().rename(spot_measure)

    G1G2_plot(Pbody_DF,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
        title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)

def G1G2_RnaProportionInPbody(Cell: pd.DataFrame, Pbody: pd.DataFrame, Spots: pd.DataFrame, spot_measure, cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          title= "Mean rna proportion in Pbodies per cell",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :

    if 'rna name' not in Spots.columns : raise MissingColumnError("'rna name' column is missing from Spots DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = update.AddCellularCycle(Spots, Cell)
    Pbody_DF = update.AddCellularCycle(Pbody, Cell)
    Spots_in_pbodies = Pbody_DF.groupby(["rna name", "cellular_cycle", "CellId"])[spot_measure].sum().rename(spot_measure)
    Spots_total = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["id"].count()
    SpotsProportion_in_pbodies = Spots_in_pbodies / Spots_total

    G1G2_plot(SpotsProportion_in_pbodies,cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
              title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)



def G1G2_CellNumber(Cell: pd.DataFrame,cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          title= "Number of cell computed", legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :

    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    cell_number = Cell.groupby(["rna name", "cellular_cycle"])["id"].count()

    G1G2_plot(cell_number,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
              title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)



def G1G2_spots_per_cell(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type: str,cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          title= "mean spot number per cell",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")
    if spots_type not in ['rna', 'malat1'] : raise ValueError("Unsupported spot type : {0}. Expected values are 'rna' or 'malat1'.".format(spots_type))

    Spots_DF = update.AddCellularCycle(Spots[Spots['spots_type'] == spots_type], Cell)
    Spots_DF = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["id"].count()
    G1G2_plot(Spots_DF,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
              title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)
    

def G1G2_cyto_spots_InPbody(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type: str, cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          title= "Mean cytoplasmic spot in Pbodies per cell",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = update.AddCellularCycle(Spots[Spots['spots_type'] == spots_type], Cell)
    Spots_DF["InCyto"] = 1 - (Spots_DF["InNucleus"].astype(bool) | Spots_DF["PbodyId"].isna())
    Spots_DF = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["InCyto"].sum()
    G1G2_plot(Spots_DF,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
             title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)



def G1G2_cyto_spots_InPbody_proportion(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type:str,cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          title= "Mean cytoplasmic spot proportion in Pbodies per cell",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = update.AddCellularCycle(Spots[Spots['spots_type'] == spots_type], Cell)
    Spots_DF["InCyto"] = 1 - ( Spots_DF["InNucleus"].astype(bool) | Spots_DF["PbodyId"].isna())
    count = Spots_DF.groupby(["rna name", "cellular_cycle", "CellId"])["InCyto"].sum()
    total = Spots_DF.query("InNucleus == 0").groupby(["rna name", "cellular_cycle", "CellId"])["id"].count()
    Spots_DF = count/total
    assert Spots_DF[Spots_DF > 1].empty, "Error : proportion > 1 found while computing cytoplasmic spots proportion in Pbody"
    G1G2_plot(Spots_DF,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
              title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)


def G1G2_total_spotnumber(Cell: pd.DataFrame, Spots : pd.DataFrame, spots_type,cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                          title= "Total spot number",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    Spots_DF = update.AddCellularCycle(Spots[Spots['spots_type'] == spots_type], Cell)
    Spots_DF = Spots_DF.groupby(['rna name', 'cellular_cycle'])['id'].count()
    G1G2_plot(Spots_DF,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
              title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)
    
def G1G2_PbodiesPerCell(Cell: pd.DataFrame, Pbody: pd.DataFrame, cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                        title= "Mean P-body number per cell",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :
    
    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    df = update.AddCellularCycle(Pbody, Cell)
    df = df.groupby(["rna name", "cellular_cycle","CellId"])["id"].count()
    G1G2_plot(df,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
              title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)
    

def G1G2_spots_per_area(Cell: pd.DataFrame, Spots: pd.DataFrame, spots_type,cellular_cycle_x = 'g1', cellular_cycle_y = 'g2',
                        title= "Mean spot number per μm²",legend=True, reset= False, close= False, show= False, path_output= None, ext ='png', **kargs) :

    if 'rna name' not in Cell.columns : raise MissingColumnError("'rna name' column is missing from Cell DF : consider using update.AddRnaName")
    if 'cellular_cycle' not in Cell.columns : raise MissingColumnError("'cellular_cycle' column is missing Cell DF : consider using update.from_IntegratedSignal_spike_compute_CellularCycleGroup")

    df = update.AddCellularCycle(Spots[Spots['spots_type'] == spots_type], Cell)
    spots_per_cell = df.groupby(["rna name", "cellular_cycle","CellId"])["id"].count()
    area_df = Cell.loc[:,['rna name', 'cellular_cycle', 'id', 'nucleus area (nm^2)']].rename(columns={'id' : 'CellId'}).set_index(['rna name', 'cellular_cycle', 'CellId']).sort_index().squeeze() / 10**6
    df = spots_per_cell / area_df
    G1G2_plot(df,  cellular_cycle_x = cellular_cycle_x, cellular_cycle_y = cellular_cycle_y,
              title= title, legend=legend, reset= reset, close= close, show= show, path_output= path_output, ext = ext, **kargs)


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