import pandas as pd
import pbwrap.plot.results_plot as rplot
import pbwrap.plot.utils as uplot
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
import CustomPandasFramework.PBody_project.update as update
from bigfish.detection import get_object_radius_nm

in_path = "/home/floricslimani/Documents/Projets/1_P_body/results_tables/"
output_path = "/home/floricslimani/Documents/Projets/1_P_body/results_tables/'"
if not in_path.endswith('/') : in_path += '/'

Acquisition = pd.read_feather(in_path + 'Acquisition')
Cell = pd.read_feather(in_path + 'Cell')
clean_Acquisition, clean_Cell = update.from_detectionthreshold_remove_acquisition(Acquisition, Cell, 80)
clean_Acquisition, clean_Cell = update.from_detectionthreshold_remove_acquisition(clean_Acquisition, clean_Cell, 250, keep='lower')

thresholds = list(Acquisition.loc[:, "RNA spot threshold"])
gene_list = list(Acquisition.loc[:, "rna name"])
# rplot.RNA_in_pbody(Acquisition, Cell, path_output= output_path + "RNA in pbody", show= False)
# rplot.RNA_in_pbody(clean_Acquisition, clean_Cell, path_output= output_path + " Clean RNA in pbody", show= False)
# rplot.RNApercentage_in_out_nucleus(Acquisition, Cell, path_output= output_path + "RNA percentage in and out Nucleus", show= False)
# rplot.RNApercentage_in_out_nucleus(clean_Acquisition, clean_Cell, path_output= output_path + "Clean RNA percentage in and out Nucleus", show= True)
# rplot.threshold(Acquisition, path_output= output_path + "Auto threshold", show= False)
# rplot.Malat_inNuc_asDapiIntensity(Cell, projtype= 'mean', path_output= output_path + "Clean Malat in fct dapi", show= False)
rplot.Malat_inNuc_asDapiIntensity(clean_Cell, projtype= 'mip', path_output= output_path + "Malat in fct dapi", show= True)