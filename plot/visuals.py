import bigfish.stack as stack
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation
from .utils import from_label_get_centeroidscoords


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



def nucleus_signal_control(Acquisitionid: int, Cell: pd.DataFrame, dapi: np.ndarray, nucleus_label: np.ndarray, measurement = "nucleus_mean_mean_signal",
                           path_output= None, show= True, axis= False, close= True):

    measures = ["id","nucleus_mean_mean_signal","nucleus_mean_median_signal","nucleus_mip_mean_signal","nucleus_mip_median_signal","plot index"]
    if measurement not in measures : raise ValueError("measurement : {0} is not in available measures : {1}".format(measurement,measures))


    centroidsdict = from_label_get_centeroidscoords(nucleus_label)
    nucleus_mask = nucleus_label.astype(bool)

    assert nucleus_label.max() == len(centroidsdict['label']), "Correct datashape should ensures that this assertion is true."


    #Figure
    fig = plt.figure(figsize=(20,20))
    implot = plt.imshow(stack.rescale(dapi))
    implot.axes.get_xaxis().set_visible(axis)
    implot.axes.get_yaxis().set_visible(axis)
    plt.tight_layout()

    #TODO This code might be update when Cell DataFrame will contain centroids coordinates.
    Cell = Cell.query("AcquisitionId == {0}".format(Acquisitionid)).loc[:, measures]
    Y = centroidsdict["centroid-0"]
    X = centroidsdict["centroid-1"]
    centroids = zip(Y,X)
    for label in centroidsdict['label'] :
        y,x = next(centroids)
        y,x = round(y), round(x)
        plt.annotate(str(label), [round(x), round(y)])
    
    if not axis : plt.cla()
    if show : plt.show()
    if path_output != None :
        stack.check_parameter(path_output = (str))
        plt.savefig(path_output)
    if close : plt.close()