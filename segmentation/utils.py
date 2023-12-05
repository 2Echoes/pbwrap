# -*- coding: utf-8 -*-
# Author: Floric Slimani <floric.slimani@live.fr>

"""
Utility functions for pbwrap.segmentation subpackage.
"""
import pandas as pd
import numpy as np
import bigfish.stack as stack
import scipy.ndimage as ndi

from skimage.measure import regionprops_table
from pbwrap.integrity import check_sameshape
from bigfish.stack import check_parameter, check_array



def distance_transform(image, label):
    """Compute distance transform of label using scipy euclidian distance transform but add weight using image pixel values.
    
    Parameters
    ----------
        image : np.ndarray
        label : np.ndarray
            Must have the same shape as image.

    Returns
    -------
    distance_transform_res : np.ndarray
    """
    stack.check_parameter(image = (np.ndarray), label = (np.ndarray))
    check_sameshape(image, label)
    distance = ndi.distance_transform_edt(label)
    distance = distance/distance.max()
    distance_transform = image * distance

    return distance_transform




def get_histogramm_highest_varation_value(array, bins= None):
    """Returns the value corresponding to the point where 1st derivative absolute value is the highest in array histogram.
    Will never return 1st elmt of the hist.

    Parameters
    ----------
        array : np.ndarray
        
    Returns
    -------
        res : int or float
    """
    stack.check_parameter(array = np.ndarray)

    if bins == None:
        if array.dtype == float :
            bins = 100 #2 decimals precision

        else:
            bins = int(array.max() - array.min())
            
    count,values = np.histogram(array, bins= bins)
    gradient = np.gradient(count)
    gradient[gradient < 0] = 0
    max_derivative_index = list(np.abs(gradient)).index(np.abs(gradient).max())
    res = values[max_derivative_index]

    return res

def auto_LoG_threshold(data: np.ndarray, bins= None) :
    """
    Looks for a threshold discriminating the spike around maximum of the LoG histogram from the highest values of the histogram.
    """
    if data.ndim != 1 : data = data.flatten()
    if bins == None:
        if data.dtype == float :
            bins = 100 #2 decimals precision

        else:
            bins = int(data.max() - data.min())



    count,values = np.histogram(data, bins= bins)
    values = values[:len(values)-1]
    spike = values[count == count.max()][0]
    gradient = np.gradient(count)
    gradient[values <= spike] = -1
    var_change = gradient[gradient > 0][0]
    change_position = list(gradient).index(var_change)
    
    return values[change_position]



def merge_channels(*channels) :
    """Merges 3D image (z,y,x) channels into a 4D image (channel,z,y,x) or 2D image channesl into a 3D image (channel,y,x).
    
    Parameters
    ----------

        *channels : np.ndarray
            3D images to merge into multi-channels image. All arrays should have the same shapes and dtypes. Should be like chan1,chan2,chan3,... . 

                
    Returns
    -------
    
        multi_channel_image : np.ndarray with shape (len(*channels), z, y, x).
            4D/3D image resulting from channels merging.
    """

    #Integrity checks
    for chan in channels : stack.check_array(chan,[2,3])
    check_sameshape(*channels)
  
    #channels merging
    dim = channels[0].ndim
    img_shape = channels[0].shape
    img_dtype = channels[0].dtype
    channel_num = len(channels)
    multi_channel_image = np.zeros(shape= [channel_num] + list(img_shape), dtype= img_dtype)
    idx_num = 0
    
    if dim == 3 :
        for chan in channels : 
            multi_channel_image[idx_num,:,:,:] = chan
            idx_num +=1

    if dim == 2 :
        for chan in channels : 
            multi_channel_image[idx_num,:,:] = chan
            idx_num +=1
    
    return(multi_channel_image)

def merge_channels_fromlists(*lists) :
    """ Merge channels from lists of 3D  or 2D images, one list corresponding to one channel.
        ch1 = [im1(z,y,x), im2(z,y,x), ... ] ch2 = [im1(z,y,x), im2(z,y,x), ... ] --> output : [im1(c,z,y,x), im2(c,z,y,x)]

    Parameters
    ----------

        *lists : List[np.ndarray]
            list of 3D/2D images to merge into multi-channels image. All arrays should have the same shapes and dtypes.

                
    Returns
    -------
    
        multi_channel_list : List[np.ndarray] 
            List of images which first axis corresponds to channels put in input.
    """
    #Integrity checks
    for lis in lists : check_parameter(lis = (list))
    check_sameshape(*lists)

    #Merging
    multi_channel_list = []
    groups = zip(*lists)
    for group in groups :
        multi_channel_list += [merge_channels(*group)]
    
    return multi_channel_list


def unstack_slices(image3D) :
    """Convert a 3D image to a list of 2D image where each images correspond to a z plane.
    
    Parameters
    ----------

        3Dimage : np.ndarray (z,y,x)
            3D image to unstack. Should be 3D with z planes being the slices to unstack.
                
    Returns
    -------
    
        slices : List[np.ndarry(y,x)].
            List of slices (z-planes) resulting from unstacking 3Dimage.
    """ 

    #Integrity checks
    check_array(image3D,3)

    #Unstacking
    slices = [slice for slice in image3D]
    return(slices)


def stack_slices(slices) :
    """Convert a list or tupple of 2D images to 3D image where each images correspond to a z plane.
    
    Parameters
    ----------

        slices : list/tuple[np.ndarray] (y,x)
            list of 2D images to stack.
                
    Returns
    -------
    
        image : np.ndarry(z,y,x).
            3Dimage resulting from slices stacking.
    """

    #Integrity
    check_parameter(slices = (list, tuple))
    for zslice in slices : check_array(zslice, ndim= 2)
    check_sameshape(*slices)
    slices = list(slices)

    #stacking

    return np.array(slices)




def euclidian_distance(pointA, pointB) :
    """Compute the euclidian distance in the plane from point A(xa ; ya) to point B(xb ; yb) : d = sqrt((xa-xb)^2 + (ya-yb)^2)
    
    Parameters
    ----------

        pointA : list[scalar]
        pointB : list[scalar]
        
    Returns
    -------
        res : float
    
    """
    #Integrity checks
    check_parameter(pointA = (list), pointB = (list))

    #Computing
    res = np.sqrt(np.square(pointA[0] - pointB[0]) + np.square(pointA[1] - pointB[1]))
    return res




def measure_Centroid(label_2D) :
    """Given a 2D labelled image, returns the coordinates (axis0, axis1) of the geometric center of each labelled regions
    
    Parameters
    ----------
        label_2D : np.ndarray(ndim = 2)
            Array containing the labeled image on which centroid measurement is performed.

    Returns
    -------
        Centroid : pd.Dataframe
            Dataframe : index = ['label', 'centroid-0', 'centroid-1']
    
        """
    #Integrity checks
    check_parameter(label_2D = (np.ndarray))


    properties_dic = regionprops_table(label_2D, properties= ["label","centroid"])
    Centroid = pd.DataFrame(properties_dic)
    return Centroid


def measure_Centroid2centroid(Centroid_previousslice, Centroid_currentslice) :
    """Measures the euclidian distance separing each centroid of {currentslice} from each centroid of {previousslice}.
    
    Parameters
    ----------
        centroid_previousslice : pd.Dataframe
            Dataframe containing the information on centroid localisation for each region. Should be computed from measure_Centroid.
        centroid_currentslice : pd.Dataframe
            Dataframe containing the information on centroid localisation for each region. Should be computed from measure_Centroid.
    
    Returns
    -------
        Centroid2centroid_measurement : pd.Dataframe
            Dataframe containing the distance and labels informations. Index = ['current slice label', 'previous slice label', 'distance (px)']
    """

    curr_label_measures = []
    prev_label_measures = []
    distance_measures = []

    #Measurement loop
    for index_currentslice in Centroid_currentslice.index :
        current_centroid = [Centroid_currentslice.at[index_currentslice, "centroid-0"], Centroid_currentslice.at[index_currentslice, "centroid-1"]]
        for index_previousslice in Centroid_previousslice.index :
            previousslice_centroid = [Centroid_previousslice.at[index_previousslice, "centroid-0"], Centroid_previousslice.at[index_previousslice, "centroid-1"]]
            distance_measures += [euclidian_distance(current_centroid, previousslice_centroid)]
            curr_label_measures += [Centroid_currentslice.at[index_currentslice, "label"]]
            prev_label_measures += [Centroid_previousslice.at[index_previousslice, "label"]]
    
    #Building output frame
    Centroid2centroid_measurement = pd.DataFrame({
        "current slice label" : curr_label_measures,
        "previous slice label" : prev_label_measures,
        "distance (px)" : distance_measures
    })
 
    return Centroid2centroid_measurement




def label_giver(Centroid_currentslice, Centroid2centroid_measurement, maximum_distance, new_label_number) :
    """Returns a data frame with current region label and new label to be assigned.
    
    Parameters
    ----------
        Centroid_currentslice : pd.Dataframe
            Dataframe containing the information on centroid localisation for each region. Should be computed from measure_Centroid.
        Centroid2centroid_measurement : pd.Dataframe
            Dataframe containing the distance and labels informations. Index = ['current slice label', 'previous slice label', 'distance (px)']
        maximum_distance : scalar
            Maximum distance between 2 centroids for label stitching.

    Returns
    -------
        label_giving : pd.Dataframe
            Index = ['current label', 'new label']

    """

    current_label = []
    new_label = []
    all_current_label = Centroid_currentslice.value_counts(subset= "label").reset_index(drop= False).drop(0, axis=1) #Group by labels

    #label giving
    Centroid2centroid_measurement = Centroid2centroid_measurement.sort_values("distance (px)").reset_index(drop= True)
    for measure in Centroid2centroid_measurement.index :
        if Centroid2centroid_measurement.at[measure, "previous slice label"] in new_label or Centroid2centroid_measurement.at[measure, "distance (px)"] > maximum_distance : continue
        current_label += [Centroid2centroid_measurement.at[measure, "current slice label"]]
        new_label += [Centroid2centroid_measurement.at[measure, "previous slice label"]]

    #building output frame
    label_giving = pd.DataFrame({
        "current label" : current_label,
        "new label" : new_label
    })

    label_giving = pd.merge(all_current_label, label_giving, how= "left", left_on= "label", right_on= "current label")
    for corres in label_giving.index :
        if not label_giving.at[corres, "new label"] >-1 : #Checking if new label is NaN
            label_giving.at[corres, "new label"] = new_label_number
            new_label_number +=1

    return label_giving, new_label_number




def relabelling(current_label, label_giving) :
    """Returns a 2D labelled image from matches between new and current labels from label_giving.
    
    Parameters
    ----------
        current_label : np.ndarray(ndim=2)
            2D labelled image which label are to be updated.
        label_giving : pd.Dataframe
            Dataframe containing matches between current and old label. Should be computed from label_giver.
            
    Returns
    -------
        new_label : np.ndarray(ndim=2)
            2D labelled image with new labels.
    
    """
    

    img_shape = current_label.shape
    new_label = np.zeros(img_shape, dtype= np.uint64)

    for region in label_giving.index :
        new_label[current_label == label_giving.at[region, "label"]] = int(label_giving.at[region, "new label"])
    
    return new_label



 

def from2Dlabel_to3Dlabel(labels, maximal_distance= 20) :
    # First of all the algorithm reconstruct the 3D image starting from z = 0. So first step is to label z=0, this label will then be used to computed z = z+1 label etc..
    # Loop : 
    #   1. Compute centroids coordinates for each label at z
    #   2. Compute centroids coordinates for each label at z-1
    #   3. Compute every euclidian distance (2D) from z = z labels to z = z-1 labels
    #   4. Sort results in increasing order.
    #   5. Going through this list label from z = z acquire label from z = z-1 :
    #       Once a label (z=z) acquires a label(z=z-1) this label cannot be acquired again
    #       A label cannot acquire a label (z=z-1) if the distance is > to 'maximal_distance' parameter in pixel
    #   6. If a label(z = z-1) is not acquired by any label (z=z) its number is deleted and will not be given during the rest of the algorithm
    #   7. If a label(z = z) finds no label(z = z -1) to acquire a new label is created.
    #   Note : If i remember correctly 6 et 7 are achieved using a 'taken_label' liste

    """
    Labels and stitches together a list of 2D mask into a 3D label uniformising labels so that object segmentation is presevered z wise.
    This operation is performed by calculating the centroid position for each layer and assigning to each region the label from the previous plane region which is closest.  
    
    Parameters
    ----------
        labels : list[np.ndarray (y,x)]
            All labelled slices must have the same shape and bool/int dtypes.
        maximal_distance : int
            A label cannot acquire a label (z=z-1) if the distance is > to 'maximal_distance' parameter in pixel
        
        
    Returns
    -------
        label3D : np.darray(z,y,x)
    
    """
    #Integrity checks
    check_parameter(labels = (list, np.ndarray), maximal_distance =(int, float))
    check_sameshape(*labels)
    import bigfish.plot as plot
    for label in labels : 
        check_array(label, ndim= 2, dtype= [np.int8, np.int16, np.uint8, np.uint16, np.int32, np.int64])

    plot.plot_images([*labels], contrast= True)

    #stitching
    
    img_shape = [len(labels)] + [labels[0].shape[0]] + [labels[0].shape[1]]

    label3D = np.zeros(img_shape, dtype= np.int64)

    label3D[0,:,:] = labels[0]
    new_label_number = int(labels[1].max()) + 1

    for z in range(1,img_shape[0]) :
        Centroid_previousslice = measure_Centroid(label3D[z-1,:,:])
        Centroid_currentslice = measure_Centroid(labels[z])
        Centroid2centroid = measure_Centroid2centroid(Centroid_previousslice, Centroid_currentslice) 
        label_giving, new_label_number = label_giver(Centroid_currentslice, Centroid2centroid, maximal_distance, new_label_number)
        label3D[z,:,:] = relabelling(labels[z], label_giving)

    return label3D