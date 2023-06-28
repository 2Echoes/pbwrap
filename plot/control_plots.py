import matplotlib.pyplot as plt
import numpy as np
import bigfish.stack as stack
import bigfish.plot as plot 
from skimage.measure import regionprops
from .utils import from_label_get_centeroidscoords, get_colors_list

"""
This submodule contains plot meant to be called during analysis pipeline.
"""

def plot_labels(labelled_image: np.ndarray, arrows = None, path_output:str = None, show= True, axis= False, close= True):
    """
    Plot a labelled image and indicate the label number at the center of each region.
    """
    #TODO : Comment
    stack.check_parameter(labelled_image = (np.ndarray, list), show = (bool))
    if isinstance(labelled_image, np.ndarray) : 
        stack.check_array(labelled_image, ndim= 2)
        labelled_image = [labelled_image]
    
    if type(arrows) == type(None) : arrows = [False] * len(labelled_image)
    else : stack.check_parameter(arrows = (list,tuple))

    plt.figure(figsize= (10,10))
    rescaled_image = stack.rescale(np.array(labelled_image[0], dtype= np.int32), channel_to_stretch= 0)
    plot = plt.imshow(rescaled_image)
    plot.axes.get_xaxis().set_visible(axis)
    plot.axes.get_yaxis().set_visible(axis)
    plt.tight_layout()

    color_list = get_colors_list(len(labelled_image))
    color_list = ['white','red']

    for index in range(0, len(labelled_image)) :
        arrow = arrows[index]
        centroid_dict = from_label_get_centeroidscoords(labelled_image[index])
        labels = centroid_dict["label"]
        Y = centroid_dict["centroid-0"]
        X = centroid_dict["centroid-1"]
        centroids = zip(Y,X)

        for label in labels :
            y,x = next(centroids)
            y,x = round(y), round(x)
            if arrow : arrow_prop = {"width" : 5, "headwidth" : 2, "headlength" : 2}
            else : arrow_prop = {}
            an = plt.annotate(str(label), [round(x), round(y)], color= color_list[index], arrowprops= arrow_prop)

    if not axis : plt.cla
    if show : plt.show()
    if path_output != None :
        stack.check_parameter(path_output = (str))
        plt.savefig(path_output)
    if close : plt.close()

    return plot

def plot_detection_steps(raw_im, spots, spots_postdecomp, cluster, contrast= True, cmap= "gray", rna_name= None, path_output= None, ext= None, show= True) :
    """Plot all detections result. Need all results from spots detection, spots decompisition and cluster detection. To plot individual detection result use bigfish.plot.
    
    Parameters
    ----------
    raw_im : np.ndarray
        A 2-d image with shape (y, x).
    spots : list or np.ndarray
        Array with coordinates and shape (nb_spots, 3) or (nb_spots, 2).
    spots_postdecom : list or np.ndarray
        Array with coordinates and shape (nb_spots, 3) or (nb_spots, 2).
    cluster : list or np.ndarray
        Array with coordinates and shape (nb_cluster, 3) or (nb_cluster, 2).
    contrast : bool
    
    """
    if contrast : raw_im= stack.rescale(raw_im, channel_to_stretch= 0)
    if not isinstance(spots, list):
        spots = [spots]
    if not isinstance(spots_postdecomp, list):
        spots_postdecomp = [spots_postdecomp]
    if not isinstance(cluster, list):
        cluster = [cluster]
    


    fig, ax = plt.subplots(2,2, sharex='col', figsize= (20,20))
    ax[0][0].imshow(raw_im, cmap= cmap)
    ax[0][1].imshow(raw_im, cmap= cmap)
    ax[1][0].imshow(raw_im, cmap= cmap)
    ax[1][1].imshow(raw_im, cmap= cmap)

    #spots
    for i, coordinates in enumerate(spots):

        # get 2-d coordinates
        if coordinates.shape[1] == 3:
            coordinates_2d = coordinates[:, 1:]
        else:
            coordinates_2d = coordinates

        # plot symbols
        y,x = zip(*coordinates_2d)
        ax[0][1].scatter(x,y, s= 1, color= 'red', linewidth= 0.5)

    for i, coordinates in enumerate(spots_postdecomp):

        # get 2-d coordinates
        if coordinates.shape[1] == 3:
            coordinates_2d = coordinates[:, 1:]
        else:
            coordinates_2d = coordinates

        # plot symbols
        y,x = zip(*coordinates_2d)
        ax[1][0].scatter(x,y, s=1, color= 'red', linewidth= 0.5)

    for i, coordinates in enumerate(cluster):

        # get 2-d coordinates
        if coordinates.shape[1] == 3:
            coordinates_2d = coordinates[:, 1:]
        else:
            coordinates_2d = coordinates

        # plot symbols
        y,x = zip(*coordinates_2d)
        ax[1][1].scatter(x,y, s=4, color= 'blue', linewidth= 1.5)


    #titles and layout
    if rna_name == None: ax[0][0].set_title("raw image",fontweight="bold", fontsize=10)
    else: ax[0][0].set_title("raw image ({0})".format(rna_name), fontweight = "bold", fontsize= 10)
    ax[0][1].set_title("spot detection : {0} spots".format(len(spots[0])),fontweight="bold", fontsize=10)
    ax[1][0].set_title("spot decomposition {0} spots".format(len(spots_postdecomp[0])),fontweight="bold", fontsize=10)
    ax[1][1].set_title("clusters detected {0} clusters".format(len(cluster[0])),fontweight="bold", fontsize=10)

    ax[0][0].axis("off")
    ax[0][1].axis("off")
    ax[1][0].axis("off")
    ax[1][1].axis("off")
    plt.tight_layout()

    # output
    if path_output is not None:
        save_plot(path_output, ext)
    if show:
        plt.show()
    else:
        plt.close()


def save_plot(path_output, ext):
    """Save the plot.

    Parameters
    ----------
    path_output : str
        Path to save the image (without extension).
    ext : str or List[str]
        Extension used to save the plot. If it is a list of strings, the plot
        will be saved several times.

    """
    # add extension at the end of the filename
    if ext == None : ext ='png'
    extension = "." + ext
    if extension not in path_output:
        path_output += extension

    # save the plot
    if isinstance(ext, str):
        # add extension at the end of the filename
        extension = "." + ext
        if extension not in path_output:
            path_output += extension
        plt.savefig(path_output, format=ext)
    elif isinstance(ext, list):
        for ext_ in ext:
            # add extension at the end of the filename
            extension = "." + ext_
            if extension not in path_output:
                path_output += extension
            plt.savefig(path_output, format=ext_)
    else:
        Warning("Plot is not saved because the extension is not valid: "
                "{0}.".format(ext))
        

def plot_cell(cell, title= None, path_output= None, show= True):
    """Plot cell plot from bigfish.plot.
    Parameters
    ----------
        cell : dict
            Result from multistack.extract cell.
        title : str
        path_output : str
        show : bool
        """
    
    cell_mask = cell["cell_mask"]
    cell_coord = cell["cell_coord"]
    nuc_mask = cell["nuc_mask"]
    nuc_coord = cell["nuc_coord"]
    rna_coord = cell["rna_coord"]
    foci_coord = cell["foci"]
    ts_coord = cell["transcription_site"]
    image_contrasted = cell["image"]

    plot.plot_cell(
            ndim=3, cell_coord=cell_coord, nuc_coord=nuc_coord, 
            rna_coord=rna_coord, foci_coord=foci_coord, other_coord= ts_coord, 
            image=image_contrasted, cell_mask=cell_mask, nuc_mask=nuc_mask, 
            title= title, show= show, path_output= path_output)
    