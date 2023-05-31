import bigfish.stack as stack
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.ndimage import binary_dilation
from skimage.segmentation import find_boundaries
from .utils import format_array_scientific_notation



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
                           use_scientific_notation= False, value_multiplicator = 1,
                           title="None", path_output= None, show= True, axis= False, close= True):


    #Figure
    fig = plt.figure(figsize=(20,20))
    implot = plt.imshow(stack.rescale(dapi), cmap= 'gray')
    implot.axes.get_xaxis().set_visible(axis)
    implot.axes.get_yaxis().set_visible(axis)
    plt.tight_layout()
    
    plot_label_boundaries(label= nucleus_label, boundary_size=boundary_size)
    if type(spots_coords) != type(None) : plot_spots(spots_coords,1)

    measures = np.array(measures, dtype= float) * value_multiplicator
    if use_scientific_notation : measures = format_array_scientific_notation(measures)
    else : measures = np.round(measures, decimals= 1)
   

    for measure, centroid in zip(measures, cells_centroids) :
        y,x = centroid
        y,x = round(y), round(x)
        plt.annotate(str(measure), [round(x), round(y)],color='black')


    plt.title(title)
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