import pandas as pd
import CustomPandasFramework.PBody_project.update as update
import bigfish.stack as stack
import pbwrap.plot.visuals as visu
import pbwrap.plot.box as box
import pbwrap.plot.scatter as scatter
import pbwrap.plot as plot

in_path = "/home/floricslimani/Documents/Projets/1_P_body/stack_O8_p21/output/20230531 17-01-21/result_tables"
output_path = "/home/floricslimani/Documents/Projets/1_P_body/Workshop/"
if not in_path.endswith('/') : in_path += '/'
if not output_path.endswith('/') : output_path += '/'

# Cleaning tables
Acquisition = pd.read_feather(in_path + 'Acquisition')
Cell = pd.read_feather(in_path + 'Cell')
rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(Acquisition, Cell, 80)
rna_clean_Acquisition, rna_clean_Cell = update.from_detectionthreshold_remove_acquisition(rna_clean_Acquisition, rna_clean_Cell, 250, keep='lower')
thresholds = list(Acquisition.loc[:, "RNA spot threshold"])
gene_list = list(Acquisition.loc[:, "rna name"])
malat_clean_Acquisition, malat_clean_Cell = update.from_malat_remove_acquisition(Acquisition, Cell, limit= 10)
new_Cell = update.from_nucleus_malat_proportion_compute_CellullarCycleGroup(malat_clean_Cell)

print(plot.get_colors_list(15))

# scatter.count_Malat_per_Cell(new_Cell, malat_clean_Acquisition)


# box.box_plot(new_Cell.loc[:, "malat1 spots in nucleus"])