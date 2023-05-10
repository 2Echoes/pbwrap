import pandas as pd
import pbwrap.plot.results_plot as rplot
import pbwrap.plot.utils as uplot
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
import CustomPandasFramework.PBody_project.update as update
from bigfish.detection import get_object_radius_nm

in_path = "/home/flo/Documents/IGH projects/1_Pbody/results/20230421_FirstAnalysis/tables/"
output_path = "/home/flo/Documents/IGH projects/1_Pbody/results/20230421_FirstAnalysis/plots/"
if not in_path.endswith('/') : in_path += '/'
if not output_path.endswith('/') : output_path += '/'

Acquisition = pd.read_feather(in_path + 'Acquisition')
Cell = pd.read_feather(in_path + 'Cell')
rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(Acquisition, Cell, 80)
rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(rna_clean_Acquisition, rna_clean_Cell, 250, keep='lower')
thresholds = list(Acquisition.loc[:, "RNA spot threshold"])
gene_list = list(Acquisition.loc[:, "rna name"])
Cell["total malat spots"] = Cell["malat1 spots in nucleus"] + Cell['malat1 spots in cytoplasm']
malat_clean_Acquisition, malat_clean_Cell = update.from_malat_remove_acquisition(Acquisition, Cell, limit= 10)



#Plots

#Threshold plot
rplot.threshold(Acquisition, path_output= output_path + "threshold plot", show= False)

#MALAT PLOTS
rplot.Malat_inNuc_asDapiIntensity(malat_clean_Cell, projtype= 'mean', path_output= output_path + "MalatSpots_DapiIntensity_plot", show= False, plot_linear_regression= True)
rplot.hist_in_nuc_malat_proportion(malat_clean_Cell, path_output = output_path + "MalatNucleus_proportion_hist", show= False, title= "Proportion of malat spots detected inside nucleus")
rplot.hist_dapi_signal(malat_clean_Cell, path_output= output_path + "Histogramme dapi signal", show= False)
rplot.hist_malat_count(malat_clean_Cell, path_output = output_path + "Malat_count_hist", show= False, title= "Malat spots count", bins= 250)

#RNA PLOTS
rplot.RNA_in_pbody(rna_clean_Acquisition, rna_clean_Cell, path_output= output_path + "mean RNA in pbody", show= False)
rplot.cytoRNA_proportion_in_pbody(rna_clean_Acquisition, rna_clean_Cell, path_output= output_path + "cytoplasmic RNA proportion in pbody", title= "cytoplasmic RNA proportion in pbody", show= False)
rplot.RNA_proportion_in_pbody(rna_clean_Acquisition, rna_clean_Cell, path_output= output_path + "RNA proportion in pbody", title= "RNA proportion in pbody", show= False)
rplot.RNApercentage_in_out_nucleus(rna_clean_Acquisition, rna_clean_Cell, path_output= output_path + "RNA percentage in and out Nucleus", show= False, plot_in_and_out_bars= False)