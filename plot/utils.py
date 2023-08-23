import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter
from math import floor, ceil
from itertools import zip_longest
import functools

def from_label_get_centeroidscoords(label: np.ndarray):
    """
    Returns
    --------
      centroid : dict{"label": list, "centroid-n": list} 
        n should be replace with 1,2.. according to the axe you wish to access."""

    check_parameter(label = (np.ndarray))
    centroid = regionprops_table(label, properties= ["label","centroid"])
    return centroid





def set_axis_ticks(axis:tuple, x_ticks_number:int = None, y_ticks_number: int= None) :
    """
    Set 'ticks_number' ticks on the plot, ticks are spaced regularly using min/max value from axis tuple.

    Parameters
    ----------
        axis : tuple (xmin, xmax, ymin, ymax)
        ticks_number : int
    """
    if not isinstance(axis, (tuple, list)) : raise TypeError("axis paremeter should be a tuple or list. It is a {0}".format(type(axis)))
    if len(axis) != 4 : raise ValueError("axis parameter should be a list containing 4 float-like : xmin, xmax, ymin, ymax.")

    xmin,xmax,ymin,ymax = axis

    #X axis
    if x_ticks_number != None :
        if xmax > 0 : last_tick = round(xmax)
        else : last_tick = 0
        if xmin < 0 : first_tick = round(xmin)
        else : first_tick = 0   
        x_ticks = np.linspace(first_tick,last_tick,x_ticks_number)
        if all(np.abs(x_ticks) > 1) : x_ticks = np.round(x_ticks)
        else : x_ticks = np.round(x_ticks, decimals= 3)
        x_ticks[0] = xmin
        x_ticks[x_ticks_number-1] = xmax
        if any(x_ticks >= 10000) : x_label = format_array_scientific_notation(x_ticks)
        elif all(x_ticks < 1) and all(x_ticks > -1) : x_label = format_array_scientific_notation(x_ticks)
        else : x_label = x_ticks
        xlocs, xlabels = plt.xticks(x_ticks,x_label)

    else : xlocs, xlabels = None,None

    #Y axis
    if y_ticks_number != None :
        last_tick = ceil(ymax) 
        y_ticks = y_ticks = np.linspace(0,last_tick,y_ticks_number)
        y_ticks[0] = floor(xmin)
        y_ticks[y_ticks_number-1] = ymax
        if any(y_ticks >= 10000) : x_label = format_array_scientific_notation(y_ticks)
        elif all(y_ticks< 1) and all(y_ticks > -1) : y_label = format_array_scientific_notation(y_ticks)
        else : y_label = y_ticks
        ylocs, ylabels = plt.xticks(y_ticks, y_label)
    else : ylocs, ylabels = None,None

    return(xlocs,xlabels,ylocs,ylabels)

    
def identity(x) :
    """
    Identity function : y = x.
    """
    return x



def save_plot(path_output, ext):
    """Save the plot.

    Parameters
    ----------
    path_output : str
        Path to save the image (without extension).
    ext : str or List[str]
        Extension used to save the plot. If it is a list of strings, the plot
        will be saved several times.
    
    Code from BigFish package.
    BSD 3-Clause License

    Copyright © 2020, Arthur Imbert
    All rights reserved.
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
        

def round_up(x) :
    x = ((x // 10)+1)*10
    return x


def round_up_bis(x,digit) :
    x = ceil(x/10**digit)*10**digit
    return x


def truncate(x, digit) :
    x = x * 10**digit
    x = ceil(x)
    x = x/10**digit
    return x
    


def format_array_scientific_notation(array: np.ndarray, precision = 2) :
    """
    Format an iterable of float into scientific notation using numpy scientific notation.
    """
    res = map(functools.partial(auto_format_float_scientific, precision= precision), array)
    return list(res)



def auto_format_float_scientific(number: float, precision: int):
    """
    Format each element from an iterable of float with more than 5 digits into scientific notation using numpy scientific notation.
    Never formats 0.

    """
    
    if number == 0 : res = 0
    elif len(str(float)) < 5 :
        res = number
    else : res = np.format_float_scientific(number,precision=precision)
    
    return res


def hist_maximum(hist:tuple) :
    highest_count = np.array(hist[0]).max()
    bins_num = list(hist[0])
    index = bins_num.index(highest_count)
    if index < len(bins_num) : index +=1
    return hist[1][index]




def get_colors_list(size:int = 62) -> np.ndarray:
    """
    Get a list of color from matplotlib.colors of length 'size'. max size 62
    """
    if not isinstance(size, int) : raise TypeError("size should be an int, it is a {0}".format(type(size)))
    if size > 62 or size < 1 : raise ValueError("Size is only supported between 1 and 62")

    red = get_red_colors()
    yellow = get_yellow_colors()
    green = get_green_colors()
    blue = get_blue_colors()
    purple = get_purple_colors()
    brown = get_brown_colors()
    pink = get_pink_colors()
    black = get_black_colors()
    grey = get_grey_colors()

    color_list = list(sum([*zip_longest(red, green, blue, black, purple, grey, yellow, brown, pink)],()))
    while None in color_list : color_list.remove(None)
    return color_list[:size]


def get_red_colors() :
    return ["#D0312D", "#990F02", "#60100B", "#7E2811", "#4E0707", "#BC544B", "#680C07"]

def get_orange_colors():
    return ["#ED7014", "#FCAE1E", "#B56727", "#BE5504", "#D67229", "#E34A27", "#FF8C00"]

def get_yellow_colors():
    return ["#D6B85A", "#DFC98A", "#C8A951", "#E7C27D", "#BDA55D", "#E4D00A", "#FFEF00"]

def get_green_colors():
    return ["#3CB043", "#3A5311", "#728C69", "#AEF359", "#5DBB63", "#028A0F", "#234F1E", "#568203", "#4CBB17", "#487800"]

def get_blue_colors():
    return ["#3944BC", "#63C5DA", "#0A1172", "#281E5D", "#1338BE", "#48AAAD", "#016064", "#2832C2", "#1F456E", "#4682B4"]

def get_purple_colors():
    return ["#A32CC4", "#7A4988", "#601A35", "#A1045A", "#663046", "#311432", "#9867C5", "#880085"]

def get_pink_colors():
    return ["#FC94AF", "#F25278", "#FE7D6A", "#FD5DA8", "#E3256B", "#FF6EC7"]

def get_brown_colors():
    return ["#4B371C", "#231709", "#795C34", "#CC7722", "#65350F", "#652A0E", "#8A3324"]

def get_black_colors():
    return ["#000000"]

def get_grey_colors():
    return ["#808080", "#373737", "#594D5B", "#3E3D53", "#9897A9", "#63645E"]




def compute_textbox_positions(position_list, text_list, character_size : int) :
    """
    """
    if len(position_list) != len(text_list) : raise ValueError("position_list and text_list arguments should have the same length.")

    textbox_position = [(xmin[0], character_size*len(text)) for xmin,text in zip(position_list,text_list)]
    return textbox_position


def is_overlapping(box1, box2) :

    xmin1,xmax1,ymin1,ymax1 = box1
    xmin2,xmax2,ymin2,ymax2 = box2

    for angle in box1 :
        x,y = angle
        if x >= xmin1 and x <= xmax1 and y >= ymin1 and y <=ymax1 :
            return True
    
    for angle in box2 :
        x,y = angle
        if x >= xmin2 and x <= xmax2 and y >= ymin2 and y <=ymax2 :
            return True
    
    return False

def compute_overlapping_list(bbox_list:list) :
    """
    Returns a list where each element is the number of bbox overlapping whith the current bbox at this position.
    """

    count_list = []

    for bbox in bbox_list :
        count = 0
        bbox_list_bis = bbox_list.copy()
        bbox_list_bis.remove(bbox)
        for bbis in bbox_list_bis :
            if is_overlapping(bbox, bbis) : count +=1
        
    count_list.append(count)
    overlapping_list = list(zip(bbox,count_list))

    return overlapping_list

def sort_overlapping_list(overlapping_list) :
    count_list, bbox_list = list(zip(*overlapping_list)) 
    df = pd.DataFrame({
        "bbox" : bbox_list,
        "overlapping_count" : count_list
    }
    )

    df.sort_values(by= "overlapping_count", ascending= False)


def random_direction() :
    """
    returns randomly -1 or 1.
    """
    roll = np.random.rand()

    if roll < 0.5 : res = -1
    else : res = 1

    return res


def random_move(bbox, length= 1) :
    xmin,xmax,ymin,ymax = bbox
    x_movement = np.random.rand()
    y_movement = 1-x_movement
    x_movement *= length
    y_movement *= length
    y_direction = random_direction()
    x_direction = random_direction()

    xmin += x_direction * x_movement
    xmax += x_direction * x_movement
    ymin += y_direction * y_movement
    ymax += y_direction * y_movement

    return (xmin,xmax,ymin,ymax)