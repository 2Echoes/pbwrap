import pandas as pd
import CustomPandasFramework.PBody_project.update as update
import bigfish.stack as stack
import pbwrap.plot.visuals as visu
import pbwrap.plot.box as box
import pbwrap.plot.scatter as scatter
import pbwrap.plot.bar as bar
import pbwrap.plot as plot
from pbwrap.quantification.CurveAnalysis import simple_linear_regression
import matplotlib.pyplot as plt
import numpy as np

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



# X = np.concatenate([np.arange(100), np.arange(100,0,-1)])
# Y = (X*np.random.randint(0,100,200)) **2 + 24
# slope,intercept = simple_linear_regression(X,Y)
# plt.plot(X,Y,'k.', label = 'random distrib')
# plt.plot(X,slope*X+intercept, 'r', label = "{0}x + {1}".format(slope,intercept))
# plt.axis('tight')
# plt.legend()
# plt.show()



scatter.DapiSignal_vs_CellNumber(Cell, integrated_signal= False)


# scatter.Malat_inNuc_asDapiIntensity(Acquisition= Acquisition.query("`rna name` in ['NF1']"), Cell= Cell, plot_linear_regression= True, title= 'NF1', s= 12)
# scatter.Malat_inNuc_asDapiIntensity(Acquisition= Acquisition.query("`rna name` in ['PABPC1']"), Cell= Cell, plot_linear_regression= True, title= 'PABPC1', s=12)



# box.box_plot(new_Cell.loc[:, "malat1 spots in nucleus"])

