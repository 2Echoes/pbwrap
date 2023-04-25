import pandas as pd
import pbwrap.plot.results_plot as rplot
import pbwrap.plot.utils as uplot
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
import CustomPandasFramework.PBody_project.update as update
from bigfish.detection import get_object_radius_nm

in_path = "/home/floricslimani/Documents/Projets/1_P_body/Whole_results/20230421_FirstAnalysis/results_tables"
output_path = "/home/floricslimani/Documents/Projets/1_P_body/Whole_results/20230421_FirstAnalysis/results_tables/plots"
if not in_path.endswith('/') : in_path += '/'
if not output_path.endswith('/') : output_path += '/'

Acquisition = pd.read_feather(in_path + 'Acquisition')
Cell = pd.read_feather(in_path + 'Cell')
clean_Acquisition, clean_Cell = update.from_detectionthreshold_remove_acquisition(Acquisition, Cell, 80)
clean_Acquisition, clean_Cell = update.from_detectionthreshold_remove_acquisition(clean_Acquisition, clean_Cell, 250, keep='lower')

thresholds = list(Acquisition.loc[:, "RNA spot threshold"])
gene_list = list(Acquisition.loc[:, "rna name"])

#FetchData
import numpy as np
clean_Cell = clean_Cell.sort_values(["AcquisitionId","id"]).reset_index(drop= True)
index = clean_Cell.sort_values('malat1 spots in nucleus').head(20).index
index = np.array(index)
select = clean_Cell.loc[index,["id","AcquisitionId","malat1 spots in nucleus"]]
select = pd.merge(select, clean_Acquisition.loc[:,["id","rna name", "rootfilename"]], how= 'left', left_on='AcquisitionId', right_on='id')
clean_Cell = pd.merge(clean_Cell, clean_Acquisition.loc[:,["id","rootfilename"]], how= 'left', left_on='AcquisitionId', right_on='id')
groupbymin = clean_Cell.loc[:,["rootfilename", "id_x"]].groupby(["rootfilename"]).min()
select = pd.merge(select, groupbymin, how= 'left', left_on= "rootfilename", right_on= "rootfilename")
select.drop(axis= 1, columns= ["id_y"])
select["plot index"] = select["id_x_x"] - select["id_x_y"]
print(select.sort_values(["rna name", "plot index"]))





#Plots
# rplot.RNA_in_pbody(Acquisition, Cell, path_output= output_path + "RNA in pbody", show= False)
# rplot.RNA_in_pbody(clean_Acquisition, clean_Cell, path_output= output_path + " Clean RNA in pbody", show= False)
# rplot.RNApercentage_in_out_nucleus(Acquisition, Cell, path_output= output_path + "RNA percentage in and out Nucleus", show= False)
# rplot.RNApercentage_in_out_nucleus(clean_Acquisition, clean_Cell, path_output= output_path + "Clean RNA percentage in and out Nucleus", show= True)
# rplot.threshold(Acquisition, path_output= output_path + "Auto threshold", show= False)
# rplot.Malat_inNuc_asDapiIntensity(Cell, projtype= 'mean', path_output= output_path + "Clean Malat in fct dapi", show= False)
rplot.Malat_inNuc_asDapiIntensity(clean_Cell, projtype= 'mip', path_output= output_path + "Malat in fct dapi", show= False)