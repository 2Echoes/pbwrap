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
Cell["total malat spots"] = Cell["malat1 spots in nucleus"] + Cell['malat1 spots in cytoplasm']
malat_clean_Acquisition, malat_clean_Cell = update.from_malat_remove_acquisition(Acquisition, Cell, limit= 10)



#Plots
# rplot.hist_in_nuc_malat_proportion(clean_Cell)
# rplot.hist_dapi_signal(Cell, path_output= output_path + "Histogramme dapi signal", show= False)
# rplot.hist_malat_count(Cell, out_nucleus= False, path_output= output_path + "Histogramme malat count in nucleus", show= False, title= "Sans filtre anti mauvais malat")
# rplot.hist_malat_count(clean_Cell, out_nucleus= False, path_output= output_path + " cleaned Histogramme malat count in nucleus", show= False, title= "Avec filtre anti mauvais malat")
# # rplot.RNA_in_pbody(Acquisition, Cell, path_output= output_path + "RNA in pbody", show= False)
# rplot.RNA_in_pbody(clean_Acquisition, clean_Cell,show= True)
# rplot.RNApercentage_in_out_nucleus(Acquisition, Cell, path_output= output_path + "RNA percentage in and out Nucleus", show= False)
# rplot.RNApercentage_in_out_nucleus(clean_Acquisition, clean_Cell, path_output= output_path + "Clean RNA percentage in and out Nucleus", show= True)
# rplot.threshold(Acquisition, path_output= output_path + "Auto threshold", show= False)
# rplot.Malat_inNuc_asDapiIntensity(Cell, projtype= 'mean', path_output= output_path + "Clean Malat in fct dapi", show= False)
rplot.Malat_inNuc_asDapiIntensity(malat_clean_Cell, projtype= 'mip', path_output= output_path + "Malat in fct dapi", show= True, plot_linear_regression= True, title = "Avec filtre malat")