import bigfish.stack as stack
import CustomPandasFramework.PBody_project.update as update
import os,re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.ndimage import binary_dilation
from skimage.segmentation import find_boundaries
from .utils import format_array_scientific_notation, save_plot



def output_spot_tiffvisual(channel, spots, path_output, dot_size = 3, rescale = True):
    """Outputs a tiff image with one channel being {channel} and the other a mask containing dots where sports are located.
    
    Parameters
    ----------
        channel : np.ndarray
            3D monochannel image
        spots :  
        path_output : str
        dot_size : int
            in pixels
    """
    stack.check_parameter(channel = (np.ndarray),spots= (list, np.ndarray), path_output = (str), dot_size = (int))
    stack.check_array(channel, ndim= [2,3])
    if channel.ndim == 3 : 
        channel = stack.maximum_projection(channel)
    if len(spots[0]) == 3 : 
        new_spots = []
        for i in range(0,len(spots)) : new_spots += [[spots[i][1], spots[i][2]]] 
        spots = new_spots

    

    spots_mask = np.zeros_like(channel)
    for spot in new_spots :
        spots_mask[spot[0], spot[1]] = 1

    
    #enlarging dots
    if dot_size > 1 : spots_mask = binary_dilation(spots_mask, iterations= dot_size-1)


    spots_mask = stack.rescale(np.array(spots_mask, dtype = channel.dtype))
    
    im = np.zeros([2] + list(channel.shape))
    im[0,:,:] = channel
    im[1,:,:] = spots_mask

    if rescale : channel = stack.rescale(channel, channel_to_stretch= 0)
    stack.save_image(im, path_output, extension= 'tif')



def nucleus_signal_control(dapi: np.ndarray, nucleus_label: np.ndarray, measures: 'list[float]' ,cells_centroids: 'list[float]',spots_coords:list = None, boundary_size = 3, 
                           use_scientific_notation= False, value_multiplicator = 1, output_spotless_copy= False,
                           title="None", path_output= None, show= True, axis= False, close= True):
    
    if path_output == None and output_spotless_copy :
        raise ValueError("Cannot output a spotless copy if no output path is given.")

    #Figure
    fig = plt.figure(figsize=(20,20))
    implot = plt.imshow(stack.rescale(dapi), cmap= 'gray')
    implot.axes.get_xaxis().set_visible(axis)
    implot.axes.get_yaxis().set_visible(axis)
    plt.tight_layout()
    plt.title(title)
    plot_label_boundaries(label= nucleus_label, boundary_size=boundary_size)
    measures = np.array(measures, dtype= float) * value_multiplicator
    if use_scientific_notation : measures = format_array_scientific_notation(measures)
    else : measures = np.round(measures, decimals= 1)
   

    for measure, centroid in zip(measures, cells_centroids) :
        y,x = centroid
        y,x = round(y), round(x)
        plt.annotate(str(measure), [round(x), round(y)],color='black')

    if type(spots_coords) != type(None) :
        if type(spots_coords) != type(None) : 
            plt.savefig(path_output + "_spotless") 
        plot_spots(spots_coords,1)



    if show : plt.show()
    if path_output != None :
        stack.check_parameter(path_output = (str))
        plt.savefig(path_output)
    if close : plt.close()


def plot_label_boundaries(label, boundary_size, color= 'blue') :
    
    #Boundaries plot
    nuc_boundaries = find_boundaries(label, mode='thick')
    nuc_boundaries = stack.dilation_filter(
        image= nuc_boundaries,
        kernel_shape= "disk",
        kernel_size= boundary_size)
    nuc_boundaries = np.ma.masked_where(
        nuc_boundaries == 0,
        nuc_boundaries)
    plt.imshow(nuc_boundaries, cmap=ListedColormap([color]))

def plot_spots(spots, color= 'red', dot_size= 1):
    
    if len(spots[0]) == 3 : 
        new_spots = []
        for i in range(0,len(spots)) : new_spots += [[spots[i][1], spots[i][2]]] 
        spots = new_spots 


    y,x = zip(*spots)
    plt.scatter(x,y, c='red', s= dot_size)



def G1_G2_labeller(result_tables_path:str, gene_list: 'list[str]', output_path:str) :
    """
    
    """


    output_path = "/home/floricslimani/Documents/Projets/1_P_body/stack_O8_p21/output/20230531 17-01-21/results_plots"
    if not result_tables_path.endswith('/') : result_tables_path += '/'
    if not output_path.endswith('/') : output_path += '/'
    os.makedirs(output_path + "G1G2visuals/", exist_ok=True)

    Acquisition = pd.read_feather(result_tables_path + 'Acquisition')
    Cell = pd.read_feather(result_tables_path + 'Cell')
    Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    
    for gene in gene_list :
        gene_Cell = Cell.query("`rna name` == '{0}'".format(gene))
    
    #Path    
    segmentation_plot_path = result_tables_path.replace("result_tables/", "steps_plots/{0}/".format(gene))
    dirlist = os.listdir(segmentation_plot_path)
    
    i = 0
    for fov in ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'] :
        print("fov : ",fov)
        acquisitionid = gene_Cell["AcquisitionId"].min() + i
        
        seg_path = None
        for file in dirlist :
            target = re.findall(".*{0}f.*{1}.*_Cell_segmentation.png".format(gene, fov), file)
            if len(target) > 0 :
                print("found : ", target)
                assert len(target) == 1, "Multiple files were found which should be impossible"
                seg_path = segmentation_plot_path + target[0]
                break
        if seg_path == None : continue
        
        _G1_G2_labelling(gene_Cell, seg_path, AcquisitionId=acquisitionid,  path_output= output_path + "G1G2_Labelling_{0}".format(fov))      
        i+=1
        print("visual saved")
    print("done")





def _G1_G2_labelling(Cell : pd.DataFrame, segmentation_plot:str, AcquisitionId:int, path_output:str) :
    """
    Add G1, G2 label to cells in the  segmentation plot.

    Parameters
    ----------
        Cell : pd.DataFrame
        segmentation_plot : str
            path to the segmentation plot on which to add labelling.
        AcquisitionId : int
            key refering to Cell["AcquisitionId"]
        path_output : str
    """
    image: np.ndarray = stack.read_image(segmentation_plot)
    df = Cell.query("`AcquisitionId` == {0}".format(AcquisitionId))

    print(image.shape)
    ax = plt.imshow(image)
    for cell, label in zip(df["cell_coordinates"], df["cellular_cycle"] ):
        ax.annotate(text = label, xy= cell)
    
    save_plot(path_output, 'png')
    plt.close()