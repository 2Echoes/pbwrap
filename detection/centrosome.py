import numpy as np
import pandas as pd
from skimage.measure import regionprops_table, regionprops


def detect_centrosome(cell, max_number_candidates = 5) :
    min_y,min_x,max_y,max_x = cell['bbox']
    cell_mask = cell["cell_mask"]
    
    #Brightest spot
    coords = []
    candidate_coords = cell["centrosome_spots"]

    if len(candidate_coords) == 0 or len(candidate_coords) > max_number_candidates : # if no spot detected or too many (likely to be noise) cell is skiped
        return np.empty(shape=(0,0), dtype=int)
    signal = cell["centrosome_mip"]
    Z,Y,X = zip(*candidate_coords)
    intensity_array = signal[Y,X]

    #1st centrosome
    max_index_1 = np.argmax(intensity_array)
    intensity_1 = intensity_array[max_index_1]
    coords.append([Z[max_index_1], Y[max_index_1], X[max_index_1]])

    #opt 2nd centrosome : if spot in a 10% range of intensity
    intensity_array[max_index_1] = 0
    max_index_2 = np.argmax(intensity_array)
    intensity_2 = intensity_array[max_index_2]
    if intensity_2 > intensity_1*0.9 :
        coords.append([Z[max_index_2], Y[max_index_2], X[max_index_2]])
    return np.array(coords, dtype= int)

    #Bigger region
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