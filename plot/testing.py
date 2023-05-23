import pandas as pd
import CustomPandasFramework.PBody_project.update as update
import bigfish.stack as stack
import pbwrap.plot.visuals as visu

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

label = stack.read_array("/home/floricslimani/Documents/Projets/1_P_body/Workshop/output/nucleus_label.npy")
dapi = stack.read_array("/home/floricslimani/Documents/Projets/1_P_body/Workshop/output/nucleus_proj.npy")

visu.nucleus_signal_control(dapi,label)