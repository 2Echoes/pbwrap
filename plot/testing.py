import pandas as pd
import pbwrap.plot.results_plot as rplot
import pbwrap.plot.utils as uplot
import matplotlib.pyplot as plt
import pbwrap.data.getdata as gdata
import CustomPandasFramework.PBody_project.update as update
from bigfish.detection import get_object_radius_nm

in_path = "/home/floricslimani/Documents/Projets/1_P_body/stack_O8_p21/output/20230418 16-17-33/result_tables/"
output_path = "/home/floricslimani/Documents/Projets/1_P_body/stack_O8_p21/output/20230418 16-17-33/"
if not in_path.endswith('/') : in_path += '/'

Acquisition = pd.read_feather(in_path + 'Acquisition')
Cell = pd.read_feather(in_path + 'Cell')


print(Acquisition)
thresholds = list(Acquisition.loc[:, "RNA spot threshold"])
gene_list = list(Acquisition.loc[:, "rna name"])
fig = rplot.RNApercentage_in_out_nlucleus(Acquisition,Cell)
plt.show()