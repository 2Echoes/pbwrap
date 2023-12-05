import bigfish.stack as stack
import bigfish.plot as plot
import bigfish.segmentation as seg
import numpy as np
from scipy import ndimage as ndi
from skimage.segmentation import watershed, random_walker
from skimage.morphology import remove_small_objects
from skimage.feature import peak_local_max
import os
from pbwrap.utils import gaussian_kernel_size, compute_anisotropy_coef
from pbwrap.segmentation.utils import get_histogramm_highest_varation_value, auto_LoG_threshold
from pbwrap.plot import histogram
from pbwrap.segmentation.custom_functions import thresholding
import napari as nap


path = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/input'
if not path.endswith('/') : path += '/'
path_out = path.replace('input','output')

voxel_size = (300,103,103) #nm to pixel scale
object_size  = (400,400,400) #nm
min_pixel_volume = 10
scale = compute_anisotropy_coef(voxel_size=voxel_size)
threshold_penalty = 1

print("files loaded from {0}".format(path))
print("files saved at {0}".format(path_out))

file_list = [path + filename for filename in os.listdir(path)]


centrosome_gaussian_size = gaussian_kernel_size(object_size, voxel_size, width= 'FWTM')
print("kernel_size = ", centrosome_gaussian_size)

for file in file_list :
    im = stack.read_image(file)
    print("shape : ", im.shape)
    im_log = stack.log_filter(im, sigma= centrosome_gaussian_size)
    flat_log = im_log.flatten()
    threshold1 = get_histogramm_highest_varation_value(flat_log)
    threshold2 = auto_LoG_threshold(flat_log) * threshold_penalty
    print('auto threhsold1 : ', threshold1)
    print('auto threhsold2 : ', threshold2)
    # histogram(flat_log.flatten())
    mask = thresholding(im_log, threshold=threshold2)
    mask_clean = remove_small_objects(mask, min_size= min_pixel_volume).astype(bool)
    im_masked = im.copy()
    im_masked[~mask_clean] = 0
    labels,_ = ndi.label(mask_clean)
    print('labels masked dtype : ', labels.dtype)



    # #Distance Watershed
    # distance = ndi.distance_transform_edt(mask_clean)
    # coords = peak_local_max(distance, footprint=np.ones((1,3,3)), labels=mask_clean)
    # peak_mask = np.zeros(distance.shape, dtype=bool)
    # peak_mask[tuple(coords.T)] = True
    # markers, _ = ndi.label(peak_mask)
    # markers[markers !=0] += 1
    # markers[im <= np.percentile(im, 80)] = 1
    # # markers[im >= np.percentile(im, 95)] = 2

    # print("starting segmentation...")
    # labels = random_walker(im, markers, beta= 130, mode= 'cg_m', spacing= scale)
    # labels[labels !=0] -= 1
    # #Intensity_watershed



    im_proj = stack.mean_projection(im)
    im_log_proj = stack.mean_projection(im_log)
    napari_viewer = nap.Viewer(title= '3D view', ndisplay= 3, show= False)
    mask_layer = napari_viewer.add_image(mask, colormap= 'red', scale= scale)
    im_layer = napari_viewer.add_image(im, opacity= 0.5, scale= scale)
    label_layer = napari_viewer.add_labels(labels, scale= scale)
    napari_viewer.show(block=True)
    # plot.plot_images([im_proj, im_log_proj, stack.maximum_projection(im_masked)], rescale= True, titles = ['raw', 'log', 'im segmented'])