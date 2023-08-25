import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter
from math import floor, ceil
from itertools import zip_longest, product
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


######## FIRST IDEA : random move until no annotations are overlapping ######


def compute_textbox_positions(position_list, text_list, character_size : int) :
    """
    """
    if len(position_list) != len(text_list) : raise ValueError("position_list and text_list arguments should have the same length.")

    textbox_position = [(xmin[0], character_size*len(text)) for xmin,text in zip(position_list,text_list)]
    return textbox_position


def collect_textbox_positions():
    pass


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

    df = df.sort_values(by= "overlapping_count", ascending= False)
    return df


def random_direction() :
    """
    returns randomly -1 or 1.
    """
    roll = np.random.rand()

    if roll < 0.5 : res = -1
    else : res = 1

    return res


def random_move(bbox, length= 1, x_direction=None, y_direction= None) :
    xmin,xmax,ymin,ymax = bbox
    x_movement = np.random.rand()
    y_movement = 1-x_movement
    x_movement *= length
    y_movement *= length

    if type(x_direction) == type(None) : x_direction = random_direction()
    if type(y_direction) == type(None) : y_direction = random_direction()

    xmin += x_direction * x_movement
    xmax += x_direction * x_movement
    ymin += y_direction * y_movement
    ymax += y_direction * y_movement

    return (xmin,xmax,ymin,ymax)

def correct_overlapping_textboxes(overlapping_list) :
    pass


############ 2nd Idea divide the plot in grid of rectangles and assign one rectangle to each annotation ###################

def compute_scale(fig, pos, text):

    ax = fig.gca()
    master_annotation = plt.annotate(text,pos)
    bbox = master_annotation.get_window_extent()
    x0,y0,x1,y1 = ax.transData.inverted().transform(bbox)

    box_xlength = x1 - x0
    box_ylength = y1 - y0

    return box_xlength, box_ylength

def compute_annotation_df(pos_list, text_list):
    annotation_df = pd.DataFrame({
        'position' : pos_list,
        'annotation' : text_list,
        'grid_coords' : np.NaN
    })

    return annotation_df

def compute_grid(x_unit, y_unit) :

    #Dividing plot area
    xmin, xmax, ymin, ymax = plt.axis()
    x_length = (xmax - xmin) // x_unit
    y_length = (ymax - ymin) // y_unit
    if (xmax - xmin) % x_unit != 0 : x_length += 1
    if (ymax - ymin) % y_unit != 0 : y_length += 1

    y_coords = np.arange(0, y_length)
    x_coords = np.arange(0, x_length)


    #Computing grid
    coordinates_list = list(product(x_coords, y_coords))
    x_coords, y_coords = zip(*coordinates_list)

    grid = pd.DataFrame({
        "coord" : coordinates_list,
        "x" : x_coords,
        "y" : y_coords,
        "empty" : [True] * len(coordinates_list)
    })

    #Excluding border
    ymax = grid["y"].max()
    xmax = grid["x"].max()
    border_idx = grid.query("y == 0 or y == {0} or x == 0 or x == {1}".format(ymax,xmax)).index
    grid.loc[border_idx, "empty"] = False

    return grid

def find_grid_coordinates_list(elmt_coords_list, x_unit, y_unit) :
    x,y = zip(*elmt_coords_list)

    x_coord = np.array(x)//x_unit
    y_coord = np.array(y)//y_unit

    return list(zip(x_coord,y_coord))

def fill_grid(coordinates_list, grid):
    x_list, y_list = zip(*coordinates_list)
    index = grid.query("x in {0} and y in {1}".format(x_list, y_list)).index
    grid.loc[index, "empty"] = False
    return grid


def find_closest_available_space(coords, grid: pd.DataFrame) :
    available_grid = grid.copy()
    taken_spaces = grid.query("empty == False").index
    available_grid = available_grid.drop(taken_spaces, axis= 0)
    x,y = coords
    available_grid["distance"] = np.sqrt( np.power((available_grid["x"] - x),2) + np.power((available_grid["y"] - x),2) )
    available_grid = available_grid.sort_values('distance').reset_index(drop= False)

    return available_grid.at[0, "index"]


def give_available_space(annotation_index, grid, annotation_df) :
    coords = annotation_df.at[annotation_index, 'position']
    space_index = find_closest_available_space(coords,grid)
    grid.at[space_index, "empty"] = False
    annotation_df.at[annotation_index, "grid_coords"] = grid.at[space_index, "coords"]

    return grid, annotation_df

def get_space_position(grid_coords, x_unit,y_unit) :

    x_pos = x_unit * grid_coords[0]
    y_pos = y_unit * grid_coords[1]

    return x_pos,y_pos



def write_annotation(annotation_df, x_unit, y_unit) :
    annotation_obj_list = []
    for idx in annotation_df.index :
        text = annotation_df.at[idx, "annotation"]
        grid_coords = annotation_df.at[idx, "grid_coords"]
        xy = get_space_position(x_unit=x_unit,y_unit=y_unit, grid_coords=grid_coords)

        annotation_obj_list.append(plt.annotate(text, xy))

    return annotation_obj_list


def annotate_plot(fig, pos_list, text_list) :
    """
    Add annotations to a plot and correct overlapping annotations.

    Parameters
    ----------
        pos_list : list[tuple]
            List of tuple, each tuple correspond to the position (x,y) of an annotation.
        text_list : list[string]
            List of string each element is the text displayed in the annotation.
    """

    if not isinstance(pos_list, list) : raise TypeError("'pos_list' argument should be a list, it is a {0}".format(type(pos_list)))
    if not isinstance(text_list, list) : raise TypeError("'text_list' argument should be a list, it is a {0}".format(type(pos_list)))
    if len(pos_list) <= 1 : raise ValueError("There should me more than 1 annotation to plot.")
    if len(pos_list) != len(text_list) : raise ValueError("pos_list and text_list should have the same number of elements.")

    x_unit, y_unit = compute_scale(fig, pos_list[0], text_list[0])
    annotation_df = compute_annotation_df(pos_list[1:],text_list[1:])
    grid = compute_grid(x_unit, y_unit)
    coords = find_grid_coordinates_list(pos_list, x_unit=x_unit, y_unit=y_unit)
    grid = fill_grid(coords, grid)

    for idx in annotation_df.index :
        grid, annotation_df = give_available_space(annotation_index= idx, annotation_df= annotation_df, grid= grid)
    annotations_obj_list = write_annotation(annotation_df,x_unit,y_unit)

    return annotations_obj_list 