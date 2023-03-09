import numpy as np
import bigfish.stack as stack
import bigfish.segmentation as seg
import cellpose.models as models
import pandas as pd
from bigfish.stack import check_array, check_parameter
from pbwrap.integrity import check_sameshape
from skimage.measure import regionprops_table 


###### Image segmentation

def Nucleus_segmentation(dapi, diameter= 150, anisotropy= 3, use_gpu= False) :
    """3D Nucleus segmentation using Cellpose from a dapi 3D grayscale image.

    Parameters
    ----------

        dapi :      np.ndarray, ndim = 3 (z,y,x). 
            The dapi should only contains the data to be analysed, prepropressing and out of focus filtering should be done prior to this operation. 
        diameter :  Int. 
            Average diameter of a nucleus in pixel. Used to rescale cellpose trained model to incoming data.
        anisotropy: Int. 
            Define the ratio between the plane (xy) resolution to the height (z) resolution. For a voxel size (300,100,100) use 3.
        use_gpu :   Bool. 
            Enable Cellpose build-in option to use GPU.
                
    Returns
    -------
    
        Nucleus_label : np.ndarray 
            With same shape as dapi. Each object has a different pixel value and 0 is background.
    """

    #Integrity checks
    check_array(dapi, ndim= 3, dtype= [np.uint8, np.uint16, np.int32, np.int64, np.float32, np.float64])
    check_parameter(diameter= (int), anisotropy= (int), use_gpu= (bool))

    #Segmentation
    nucleus_model = models.Cellpose(gpu= use_gpu, model_type = "nuclei")
    channels = [0,0]
    nucleus_label = nucleus_model.eval(dapi, diameter= diameter, channels = channels, anisotropy= anisotropy, do_3D= True, stitch_threshold = 0.1)[0]
    for label_num in range(0,len(nucleus_label)) :
        min_objct_size = int(round((np.pi * (diameter/2)**2) /2)) # area in pixel
        nucleus_label[label_num] = seg.clean_segmentation(nucleus_label[label_num], small_object_size= min_objct_size)
    
    nucleus_label = seg.remove_disjoint(nucleus_label)


    return nucleus_label




def Cytoplasm_segmentation(cy3, dapi= None, diameter= 250, maximal_distance= 100, use_gpu= False) :
    """Due to bad performance using 3D cellpose with cy3 channel image. A 2D cell segmentation is performed for each slice and a 3D labelling is performed using a closest centroid method.

    Parameters
    ----------

        cy3 :      np.ndarray, ndim = 3 (z,y,x). 
            cy3 should only contains the data to be analysed, prepropressing and out of focus filtering should be done prior to this operation. 
        diameter :  Int. 
            Average diameter of a cell in pixel. Used to rescale cellpose trained model to incoming data.
        use_gpu :   Bool. 
            Enable Cellpose build-in option to use GPU.
    
    Returns
    -------
    
        cytoplasm_labels : List[np.ndarray] 
            List of numpy arrays with shape(x,y). Eache element correspond to a z plane and each object has a different pixel value and 0 is background.
    """
    #Integrity checks
    check_array(cy3, ndim= [2,3], dtype= [np.uint8, np.uint16, np.int32, np.int64, np.float32, np.float64])
    check_parameter(diameter= (int), dapi= (np.ndarray, type(None)), use_gpu= (bool))

    #Segmentation
    if cy3.ndim == 3 : cytoplasm_slices = unstack_slices(cy3)
    cytoplasm_model = models.Cellpose(gpu= use_gpu, model_type = "cyto")
 

    if type(dapi) == type(None) : 
        channels = [0,0]
        if cy3.ndim == 3 : image = cytoplasm_slices 
        else : image = cy3
    else :
        check_sameshape(cy3, dapi)
        channels = [1,2]
        if dapi.ndim == 3 :
            nucleus_slices = unstack_slices(dapi)
            image = merge_channels_fromlists(cytoplasm_slices, nucleus_slices)
        else :
            image = merge_channels(cytoplasm_slices,dapi)


    cytoplasm_labels = cytoplasm_model.eval(image, diameter= diameter, channels= channels, do_3D= False)[0]
    if cy3.ndim == 3 : cytoplasm_label = from2Dlabel_to3Dlabel(cytoplasm_labels, maximal_distance= maximal_distance)
    else : cytoplasm_label = cytoplasm_labels


    return cytoplasm_label




def pbody_segmentation(egfp, sigma = 2, threshold = 180, small_obj_sz= 250, fill= True) :
    """Performs Pbody segmentation on 2D or 3D egfp numpy array
    
    Parameters
    ----------
        egfp : np.ndarray(y,x) or (z,y,x)
        sigma : scalar
            Parameter applied to Gaussian filter
        thresold : int
            Parameter applied to thresholding
            
    Returns
    -------
        Pbody_label : np.ndarray(y,x) or (z,y,x)
    """

    check_parameter(egfp = (np.ndarray), sigma = (int, float), threshold = (int))
    check_array(egfp, ndim= [2,3])
    dim = egfp.ndim


    egfp_log = stack.log_filter(egfp,sigma)

    if dim == 2 :
        egfp_label = _pbody_slice_seg(egfp_log, threshold=threshold, small_object_size=small_obj_sz, fill_holes=fill)

    else :
        slices = unstack_slices(egfp_log)
        label_list = []
        for z in range(0,len(slices)):
            label_list += [_pbody_slice_seg(slices[z], threshold=threshold, small_object_size=small_obj_sz, fill_holes=fill)]
        egfp_label = from2Dlabel_to3Dlabel(label_list, maximal_distance= 50)

    return egfp_label



def _pbody_slice_seg(image, threshold, small_object_size, fill_holes) : 
        
        egfp_mask = seg.thresholding(image, threshold)
        egfp_mask = seg.clean_segmentation(egfp_mask,small_object_size= small_object_size, fill_holes= fill_holes )
        egfp_label = seg.label_instances(egfp_mask)

        return egfp_label





####### Image Operations
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
    slices = []
    for slice in image3D : slices += [slice]
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
    img_shape = [len(slices)] + list(slices[0].shape)
    img_dtype = slices[0].dtype
    image = np.zeros(img_shape, dtype= img_dtype)
    z_index = 0
    for zslice in slices:
        image[z_index,:,:] = zslice
        z_index +=1

    return image




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
    """
    Labels and stitches together a list of 2D mask into a 3D label uniformising labels so that object segmentation is presevered z wise.
    This operation is performed by calculating the centroid position for each layer and assigning to each region the label from the previous plane region which is closest.  
    
    Parameters
    ----------
        mask : list[np.ndarray (y,x)]
            All mask must have the same shape and bool/int dtypes.
        
        
    Returns
    -------
        label3D : np.darray(z,y,x)
    
    """
    #Integrity checks
    check_parameter(labels = (list), maximal_distance =(int, float))
    check_sameshape(*labels)
    for label in labels : check_array(label, ndim= 2, dtype= [np.int8, np.int16, np.uint8, np.uint16, np.int32, np.int64])

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