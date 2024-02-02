import bigfish.stack as stack
import CustomPandasFramework.PBody_project.update as update
import os,re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pbwrap.data as data
import pbwrap.segmentation as segmentation
import tifffile as tif
from matplotlib.colors import ListedColormap
from scipy.ndimage import binary_dilation
from skimage.segmentation import find_boundaries
from scipy.ndimage import distance_transform_edt
from .utils import format_array_scientific_notation, save_plot
import warnings


def output_spot_tiffvisual(channel,spots_list, path_output, dot_size = 3, rescale = True):
    
    """
    Outputs a tiff image with one channel being {channel} and the other a mask containing dots where sports are located.
    
    Parameters
    ----------
        channel : np.ndarray
            3D monochannel image
        spots : list[np.ndarray] or np.ndarray
            Spots arrays are ndarray where each element corresponds is a tuple(z,y,x) corresponding to 3D coordinate of a spot
            To plot different spots on different channels a list of spots ndarray can be passed. 
        path_output : str
        dot_size : int
            in pixels
    """
    
    stack.check_parameter(channel = (np.ndarray), spots_list= (list, np.ndarray), path_output = (str), dot_size = (int))
    stack.check_array(channel, ndim= [2,3])
    if isinstance(spots_list, np.ndarray) : spots_list = [spots_list]

    if channel.ndim == 3 : 
        channel = stack.maximum_projection(channel)

    im = np.zeros([1 + len(spots_list)] + list(channel.shape))
    im[0,:,:] = channel

    for level in range(len(spots_list)) :
        if len(spots_list[level]) == 0 : continue
        else :
            spots_mask = np.zeros_like(channel)
            
            #Unpacking spots
            if len(spots_list[level][0]) == 2 :
                Y,X = zip(*spots_list[level])
            elif len(spots_list[level][0]) == 3 :
                Z,Y,X = zip(*spots_list[level])
                del Z
            else :
                Z,Y,X,*_ = zip(*spots_list[level])
                del Z,_
            
            #Reconstructing signal
            spots_mask[Y,X] = 1
            if dot_size > 1 : spots_mask = binary_dilation(spots_mask, iterations= dot_size-1)
            spots_mask = stack.rescale(np.array(spots_mask, dtype = channel.dtype))
            im[level + 1] = spots_mask

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



def G1_G2_labeller(result_tables_path:str, grouping, input_path:str, output_path:str, gene_list:'list[str]'=None,**function_kargs) :
    """
    
    """

    if not result_tables_path.endswith('/') : result_tables_path += '/'
    if not input_path.endswith('/') : input_path += '/'
    if not output_path.endswith('/') : output_path += '/'
    output_path += "G1G2visuals/"
    os.makedirs(output_path , exist_ok=True)

    Acquisition = pd.read_feather(result_tables_path + 'Acquisition')
    Cell = pd.read_feather(result_tables_path + 'Cell')
    Cell = update.JoinCellAcquisition(Acquisition, Cell, Acquisition_columns= ["rna name"])
    
    if type(gene_list) == type(None) : gene_list = data.from_Acquisition_get_rna(Acquisition)
    print(len(gene_list), " genes found.")
    
    for gene in gene_list :
        path = output_path + "{0}/".format(gene)
        os.makedirs(path, exist_ok= True)
        gene_Cell_index = Cell.query("`rna name` == '{0}'".format(gene)).index
        gene_Cell = grouping(Cell.loc[gene_Cell_index,:], **function_kargs)

    
    #Path    
        segmentation_plot_path = result_tables_path.replace("result_tables/", "steps_plots/{0}/".format(gene))
        dirlist = os.listdir(segmentation_plot_path)

        i = 0
        for fov in ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'] :
            print("fov : ",fov)
            acquisitionid = gene_Cell["AcquisitionId"].min() + i

            seg_path = None
            for file in dirlist :
                target = re.findall("(.*{0}f.*{1}.*)_Cell_segmentation.png".format(gene, fov), file)
                if len(target) > 0 :
                    print("found : ", target)
                    assert len(target) == 1, "Multiple files were found which should be impossible"
                    print("initial target : ",target)
                    target = target[0].replace("--","-DAPI-") + ".tiff"
                    print("corrected target : ", target)
                    seg_path = input_path + target
                    break
            if seg_path == None :
                if acquisitionid != gene_Cell["AcquisitionId"].min() : i+=1 
                continue

            _G1_G2_labelling(gene_Cell, seg_path, AcquisitionId=acquisitionid,  path_output= path + "{1}_G1G2_Labelling_{0}".format(fov,gene))
            i+=1
            print("visual saved")
    print("done")


def merge_to_RGB(red, green, blue, rescale= True) :

    if red.shape != green.shape or red.shape != blue.shape :
        raise ValueError("All images to merge should have the same shape.")

    shape = red.shape

    red_im = stack.cast_img_uint8(red)
    green_im = stack.cast_img_uint8(green)
    blue_im = stack.cast_img_uint8(blue)

    if rescale :
        red_im = stack.rescale(red_im)
        blue_im = stack.rescale(blue_im)
        green_im = stack.rescale(green_im)
    
    red_im = red_im.flatten()
    blue_im = blue_im.flatten()
    green_im = green_im.flatten()

    image_rgb = zip(red_im, green_im, blue_im)
    image_rgb = np.array(list(image_rgb)).reshape(shape[0],shape[1],3)

    return image_rgb


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
    image_DAPI: np.ndarray = stack.read_image(segmentation_plot)
    image_Cy3 = stack.read_image(segmentation_plot.replace("DAPI", "Cy3"))
    image_EGFP = stack.read_image(segmentation_plot.replace("DAPI", "EGFP"))
    df = Cell.query("`AcquisitionId` == {0}".format(AcquisitionId))
    
    if image_DAPI.ndim == 3 : 
        im_shape = image_DAPI[0].shape
        image_EGFP = stack.rescale(stack.cast_img_uint8(stack.maximum_projection(image_EGFP)), channel_to_stretch= 0).flatten()
        image_Cy3 = stack.rescale(stack.cast_img_uint8(stack.maximum_projection(image_Cy3))).flatten()
        image_DAPI = stack.rescale(stack.cast_img_uint8(stack.mean_projection(image_DAPI))).flatten()
    else : im_shape = image_DAPI.shape

    image_rgb = zip(image_EGFP, image_Cy3, image_DAPI)
    image_rgb = np.array(list(image_rgb)).reshape(im_shape[0],im_shape[1],3)

    fig = plt.figure(figsize= (round(im_shape[0]/100), round(im_shape[1]/100)))
    ax = plt.imshow(image_rgb)
    plt.axis(False)
    fig.tight_layout()
    for cell, label in zip(df["cell_coordinates"], df["cellular_cycle"] ):
        plt.annotate(text = label, xy= (cell[1],cell[0]), color= 'red', size= 'large')
    save_plot(path_output, 'png')
    plt.close()


def reconstruct_boolean_signal(image_shape, spot_list: list) :
    signal = np.zeros(image_shape, dtype= bool)
    
    if len(spot_list) > 0 :
        Z, Y, X = list(zip(*spot_list))
    else : 
        return signal
    
    if len(image_shape) == 2 :
        signal[Y,X] = 1
    else :
        signal[Z,Y,X] = 1

    return signal

def dapi_artifact_overlay(dapi_channel, spots_array, path_out):
    
    if len(spots_array) != 0 :
        proj = stack.maximum_projection(dapi_channel)
        nucleus_mask = segmentation.Nucleus_segmentation(proj, use_gpu= True, model_type= 'nuclei').astype(bool)
        artifact_signal = proj.copy()
        nucleus_signal = np.zeros_like(proj)
        artifact_signal[nucleus_mask] = 0
        nucleus_signal[nucleus_mask] = 1
        artifact_signal = stack.rescale(artifact_signal, channel_to_stretch= 0)
        nucleus_signal = stack.rescale(nucleus_signal, channel_to_stretch= 0)
        
        RGB = np.zeros(dtype = dapi_channel.dtype, shape= (*proj.shape, 3))
        Z, Y, X = list(zip(*spots_array))
        RGB[:,:,0][Y, X] = 1
        RGB[:,:,0] = stack.rescale(RGB[:,:,0])
        RGB[:,:,1] = stack.rescale(artifact_signal)
        RGB[:,:,2] = stack.rescale(nucleus_signal)
        stack.save_image(RGB, path=path_out, extension= 'tiff')

def colocalisation_plot(background, shape, voxel_size, colocalisation_distance, path_output, spot_list1, spot_list2, spot_list3=[], dot_size= 2) :
    """
    Only works for 3d spots.
    """

    signal1 = reconstruct_boolean_signal(shape, spot_list1)
    signal2 = reconstruct_boolean_signal(shape, spot_list2)
    signal3 = reconstruct_boolean_signal(shape, spot_list3)
    coord_grid = np.indices(shape)

    #Coloc 1 & 2 and 1 & 3
    distance_map_signal_1 = distance_transform_edt(
        input= np.logical_not(signal1),
        sampling= voxel_size
        )
    
    distance_map_signal_1 = distance_map_signal_1 < colocalisation_distance

    mask_coloc_1_2 = np.logical_and(signal2, distance_map_signal_1)
    index_coloc_1_2 = (coord_grid[0][mask_coloc_1_2], coord_grid[1][mask_coloc_1_2], coord_grid[2][mask_coloc_1_2])
    del mask_coloc_1_2

    mask_coloc_1_3 = np.logical_and(signal3, distance_map_signal_1)
    index_coloc_1_3 = (coord_grid[0][mask_coloc_1_3], coord_grid[1][mask_coloc_1_3], coord_grid[2][mask_coloc_1_3])
    del mask_coloc_1_3

    del distance_map_signal_1

    #Coloc 2 & 3
    distance_map_signal_2 = distance_transform_edt(
        input= np.logical_not(signal2),
        sampling= voxel_size
        )
    
    distance_map_signal_2 = distance_map_signal_2 < colocalisation_distance

    mask_coloc_2_3 = np.logical_and(signal3, distance_map_signal_2)
    index_coloc_2_3 = (coord_grid[0][mask_coloc_2_3], coord_grid[1][mask_coloc_2_3], coord_grid[2][mask_coloc_2_3])
    del mask_coloc_2_3, coord_grid

    spots_colloc_1_2 = list(zip(*index_coloc_1_2))
    spots_colloc_1_3 = list(zip(*index_coloc_1_3))
    spots_colloc_2_3 = list(zip(*index_coloc_2_3))

    colloc_1_2_signal = reconstruct_boolean_signal(shape, spots_colloc_1_2)
    colloc_1_3_signal = reconstruct_boolean_signal(shape, spots_colloc_1_3)
    colloc_2_3_signal = reconstruct_boolean_signal(shape, spots_colloc_2_3)

    image = np.zeros((1,) + shape[:1] +(7,) + shape[1:])
    image[0,:,0,:,:] = background
    for channel, signal in enumerate([signal1, signal2, signal3, colloc_1_2_signal, colloc_1_3_signal, colloc_2_3_signal]) :
        image[0,:,channel+1,:,:] = binary_dilation(signal, iterations= dot_size-1)

    image = image.astype(np.int16)

    warnings.simplefilter("ignore")
    try :
        tif.imwrite(
            path_output,
            image,
            imagej= True,
            bigtiff= True,
            metadata={
             'axes': 'TZCYX',
            }
            )
    except Exception as error : raise error
    finally : warnings.simplefilter("default")
    