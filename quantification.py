import bigfish.stack as stack
import bigfish.multistack as multistack
import bigfish.plot as plot


## Testing
image_input = "/Users/floric.slimani/Documents/Projets/1_PBodies/RNA-grouped_3D-z-stacked/input/"
output = "/Users/floric.slimani/Documents/Projets/1_PBodies/RNA-grouped_3D-z-stacked/output/quantification/"
segmentation_input = "/Users/floric.slimani/Documents/Projets/1_PBodies/RNA-grouped_3D-z-stacked/seg_res/"

ndim = 3
cell_label = stack.read_array(segmentation_input + "cytoplasm_2Dlabel.npy")
nuc_label = stack.read_array(segmentation_input + "nucleus_3Dlabel.npy")
rna_coords = stack.read_array(segmentation_input + "rnaspots.npy")
malat1_coords = stack.read_array(segmentation_input + "malat1spots.npy")
other_coords = {"malat1" : malat1_coords}

print("cell label", cell_label.dtype, cell_label.shape)
print("nuc label", nuc_label.dtype, nuc_label.shape)
print("rna", rna_coords.dtype, rna_coords.shape)
print("malat1", malat1_coords.dtype, malat1_coords.shape)




fov_results = multistack.extract_cell(cell_label=cell_label[5], ndim=ndim, nuc_label= nuc_label[5], rna_coord= rna_coords, others_coord= other_coords)
for i, cell_results in enumerate(fov_results):
    print("cell {0}".format(i))
    
    # get cell results
    cell_mask = cell_results["cell_mask"]
    cell_coord = cell_results["cell_coord"]
    nuc_mask = cell_results["nuc_mask"]
    nuc_coord = cell_results["nuc_coord"]
    rna_coord = cell_results["rna_coord"]
    #foci_coord = cell_results["foci"]
    #ts_coord = cell_results["transcription_site"]
    #image_contrasted = cell_results["image"]
    print("\r number of rna {0}".format(len(rna_coord)))
    #print("\r number of foci {0}".format(len(foci_coord)))
    #print("\r number of transcription sites {0}".format(len(ts_coord)))
    
    # plot cell
    plot.plot_cell(
        ndim=3, cell_coord=cell_coord, nuc_coord=nuc_coord, 
        rna_coord=rna_coord, cell_mask=cell_mask, nuc_mask=nuc_mask, 
        title="Cell {0}".format(i))

