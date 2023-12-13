import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, regionprops


def detect_centrosome(cell, centrosome_presegmentation, size_threshold= None) :
    min_y,min_x,max_y,max_x = cell['bbox']
    cell_mask = cell["cell_mask"]
    centrosome_label_crop = centrosome_presegmentation.copy()[:, min_y: max_y, min_x : max_x]
    if cell_mask.ndim == 2 : cell_mask = np.repeat(np.expand_dims(cell_mask, axis=0),
                                                  len(centrosome_presegmentation),
                                                  axis= 0)
    centrosome_label_crop[~cell_mask] = 0
    regions_df = pd.DataFrame(regionprops_table(centrosome_label_crop, properties= ['label', 'centroid', 'area_filled']))
    regions_df = regions_df.sort_values('area_filled', ascending= False).reset_index(drop= True)


    if len(regions_df) != 0 : 
        coords = [(round(regions_df.at[0, "centroid-0"]), round(regions_df.at[0, "centroid-1"]), round(regions_df.at[0, "centroid-2"]))]
    else : return np.empty(shape=(0,0), dtype=int)
    if len(regions_df) >1 :
        if type(size_threshold) == type(None) :
            size_threshold =  max(regions_df.at[0,'area_filled'] / 2, 10)
        if regions_df.at[1, "area_filled"] >= size_threshold : 
            coords.append((round(regions_df.at[1, "centroid-0"]), round(regions_df.at[1, "centroid-1"]), round(regions_df.at[1, "centroid-2"])))

    return np.array(coords, dtype= int)



# #testing
# import os
# import bigfish.stack as stack
# import bigfish.plot as plot
# import pbwrap.segmentation as seg
# from pbwrap.utils import gaussian_kernel_size, compute_anisotropy_coef

# path = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/input'
# if not path.endswith('/') : path += '/'
# path_out = path.replace('input','output')
# im_path = path + os.listdir(path)[0]
# voxel_size = (300,103,103) #nm to pixel scale
# object_size  = (400,400,400) #nm
# min_pixel_volume = 10
# scale = compute_anisotropy_coef(voxel_size=voxel_size)
# threshold_penalty = 1


# im = stack.read_image(im_path)
# prelabel = seg.centrosome_segmentation_candidate_regions(
#     image= im,
#     centrosome_size= object_size,
#     voxel_size=voxel_size
# )


# cell = {
#     'bbox' : [593, 734, 593 + 268, 734 +374],
#     'nuc_mask' : np.ones((15, 268,374), dtype=bool)
# }

# detect_centrosome(cell=cell, centrosome_presegmentation= prelabel)