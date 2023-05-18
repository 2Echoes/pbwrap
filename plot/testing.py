import pandas as pd
import pbwrap.plot.results_plot as rplot
import pbwrap.plot.utils as uplot
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
import CustomPandasFramework.PBody_project.update as update
from bigfish.detection import get_object_radius_nm

in_path = "/home/floricslimani/Documents/Projets/1_P_body/Whole_results/20230421_FirstAnalysis/results_tables"
output_path = "/home/floricslimani/Documents/Projets/1_P_body/Whole_results/20230421_FirstAnalysis/results_tables"
if not in_path.endswith('/') : in_path += '/'
if not output_path.endswith('/') : output_path += '/'

# Cleaning tables
Acquisition = pd.read_feather(in_path + 'Acquisition')
Cell = pd.read_feather(in_path + 'Cell')
rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(Acquisition, Cell, 80)
rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(rna_clean_Acquisition, rna_clean_Cell, 250, keep='lower')
thresholds = list(Acquisition.loc[:, "RNA spot threshold"])
gene_list = list(Acquisition.loc[:, "rna name"])
Cell["total malat spots"] = Cell["malat1 spots in nucleus"] + Cell['malat1 spots in cytoplasm']
malat_clean_Acquisition, malat_clean_Cell = update.from_malat_remove_acquisition(Acquisition, Cell, limit= 10)

#Test violin plots
# rplot.violin_rna_in_pbody(Acquisition, Cell)










#Testing cytoRNAproportion
# join_frame = pd.merge(Cell, Acquisition.loc[:,["id", "rna name"]], how= "left", left_on= "AcquisitionId", right_on= "id")
# join_frame = join_frame.drop(axis= 0, index= join_frame[join_frame["pbody number"] == 0].index)
# gene_Cell = join_frame [join_frame["rna name"] == "ZNF100"]
# gene_Cell.loc[:,["cyto RNA proportion in p-body"]] = gene_Cell.loc[:, "rna spots in body"] / gene_Cell.loc[:, "nb_rna_out_nuc"]
# print(gene_Cell.loc[:,["rna spots in body","nb_rna_out_nuc", "cyto RNA proportion in p-body", "proportion_rna_centrosome", ]])