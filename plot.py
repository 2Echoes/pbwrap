import matplotlib.pyplot as plt
import numpy as np
import bigfish.stack as stack
from skimage.measure import regionprops

def plot_labels(labelled_image, path_output= None, show= True, axis= False, close= True):
    plt.figure(figsize= (10,10))
    stack.check_parameter(labelled_image = (np.ndarray), show = (bool))
    rescaled_image = stack.rescale(labelled_image, channel_to_stretch= 0)
    plot = plt.imshow(rescaled_image)
    plot.axes.get_xaxis().set_visible(axis)
    plot.axes.get_yaxis().set_visible(axis)
    plt.tight_layout()

    regions = regionprops(label_image=labelled_image)
    for props in regions :
        an = plt.annotate(str(labelled_image[round(props.centroid[0]),round(props.centroid[1])]), [round(props.centroid[1]), round(props.centroid[0])])

    if not axis : plt.cla
    if show : plt.show()
    if path_output != None :
        stack.check_parameter(path_output = (str))
        plt.savefig(path_output)
    if close : plt.close()

    return plot