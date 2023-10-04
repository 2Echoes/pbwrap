import bigfish.stack as stack
import CustomPandasFramework.PBody_project.update as update
import os,re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pbwrap.data as data
import pbwrap.segmentation as segmentation
from matplotlib.colors import ListedColormap
from scipy.ndimage import binary_dilation
from skimage.segmentation import find_boundaries
from .utils import format_array_scientific_notation, save_plot


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


def output_spot_tiffvisual_old(channel, spots, path_output, dot_size = 3, rescale = True):
    #Obselete delete if pipeline is working
    """
    Outputs a tiff image with one channel being {channel} and the other a mask containing dots where sports are located.
    
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
        for i in range(0, len(spots)) : new_spots += [[spots[i][1], spots[i][2]]] 
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
    Z, Y, X = list(zip(*spot_list))
    if len(image_shape) == 2 :
        signal[Y,X] = 1
    else :
        signal[Z,Y,X] = 1

    return signal


def fullcheck_visual(result_tables_path:str, grouping, input_path:str, output_path:str, gene_list:'list[str]'=None, pbody_kernel_size= 2.25) :
    
    """
    Output
    ------
        Tiff format output, saved at `output_path`. Image is multichannel and should contain :
            DAPI channel : with g1/g2 identification
            Cy3 channel : raw channel
            Spots detected channel
            EGFP channel
            Pbody_mask

    """
    # Loading results from pipeline
    Acquisition = pd.read_feather(result_tables_path + "Acquisition")
    Cell = pd.read_feather(result_tables_path + "Cell")
    Spots = pd.read_feather(result_tables_path + "Spots")
    
    if type(gene_list) == type(None) : gene_list = data.from_Acquisition_get_rna(Acquisition)
    print(len(gene_list), " genes found.")
    for gene in gene_list :
        path = output_path + "FullCheck/{0}/".format(gene)
        os.makedirs(path, exist_ok= True)
        Acquisition_index = Acquisition.query("`rna name` == '{0}'".format(gene)).index
        gene_Cell = pd.merge(Acquisition.loc[Acquisition_index, "id"], Cell, how= 'inner', left_on='id', right_on= 'AcquisitionId')
        gene_Cell = grouping(gene_Cell)
        for acquisition in Acquisition_index :
            segmentation_plot_path = result_tables_path.replace("result_tables/", "steps_plots/{0}/".format(gene))
            dirlist = os.listdir(segmentation_plot_path)

            for fov in ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'] :
                print("fov : ",fov)

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
                if seg_path != None : break
            #loading images
            dapi = stack.mean_projection(stack.read_image(seg_path))
            egfp = stack.mean_projection(stack.read_image(seg_path.replace('DAPI', 'EGFP')))
            cy3 = stack.maximum_projection(stack.read_image(seg_path.replace('DAPI', 'Cy3')))

            shape = dapi.shape

            egfp_rescale = stack.rescale(egfp)
            pbody_mask = segmentation.pbody_segmentation(egfp_rescale, sigma= pbody_kernel_size , threshold= None, small_object_size= 10, fill_holes= True).astype(bool)
            del egfp_rescale

            spots = list(Spots.query("spots_type == 'rna' and AcquisitionId == {0}".format(Acquisition.at[acquisition,"id"]))["spots_coords"])
            spots_signal = reconstruct_boolean_signal(shape, spots)

            # #multichannl :
            # im_zip = (
            # dapi.flatten(),
            # cy3.flatten(),
            # spots_signal.flatten(),
            # egfp.flatten(),
            # pbody_mask.flatten(),
            # )
            # image = np.array(list(im_zip)).reshape(1,shape[0],shape[1],5)
            # image = np.stack([dapi, cy3, spots_signal, egfp, pbody_mask]).reshape(5,1,shape[0],shape[1])
            dtype = dapi.dtype
            print(dtype)
            image = np.zeros(shape= [5,shape[0],shape[1]], dtype= dtype)
            image[0] = dapi
            image[1] = cy3.astype(dtype)
            # image[2,0] = spots_signal.astype(dtype)
            # image[3,0] = egfp.astype(dtype)
            # image[4,0] = pbody_mask.astype(dtype)
            # print(image.shape)
            stack.save_image(image, path= path + "Check_Visuals_{0}".format(gene + 'f' + fov), extension= 'tiff')
            im = stack.read_image(path + "Check_Visuals_{0}.tiff".format(gene + 'f' + fov))
            print(im.shape)
            quit()
            print('fov saved')