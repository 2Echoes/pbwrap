import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter
from math import floor, ceil
from itertools import zip_longest, product
import functools, itertools

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

def get_markers_generator() :
    markers = (marker for marker in ['o','v','^','<','>','s','8','p','*','h','P','D'])
    gen = itertools.cycle(markers)
    return gen

def get_hatch_generator():
    hatchs = ['/','|','-','+','*']
    gen = itertools.cycle((marker for marker in hatchs))
    return gen

def get_colors_list(size:int = 100) -> list:
    """
    Get a list of color from matplotlib.colors of length 'size'.
    100 different shade in the library
    """
    if not isinstance(size, int) : raise TypeError("size should be an int, it is a {0}".format(type(size)))
    if size < 1 : raise ValueError("Size should be >= 1")

    red = _get_red_colors()
    yellow = _get_yellow_colors()
    green = _get_green_colors()
    blue = _get_blue_colors()
    purple = _get_purple_colors()
    brown = _get_brown_colors()
    pink = _get_pink_colors()
    orange = _get_orange_colors()
    black = _get_black_colors()
    grey = _get_grey_colors()

    color_list = list(sum([*zip_longest(red, green, blue, black, orange, purple, grey, yellow, brown, pink)],()))
    length = len(color_list)
    while None in color_list : color_list.remove(None)
    if size > length :
        iteration = ceil(size / length)
        return (color_list * iteration)[:size]
        # print("length :" , length)
        # print("len returned : ", len(color_list * iteration))
        # print("size :" , size)
        # print("iteration :", iteration)return (color_list * iteration)[:size]

    return (color_list)[:size]


def _get_red_colors() :
    return ["#D0312D", "#990F02", "#60100B", "#7E2811", "#4E0707", "#BC544B", "#680C07"]

def _get_orange_colors():
    return ["#ED7014", "#FCAE1E", "#B56727", "#BE5504", "#D67229", "#E34A27", "#FF8C00"]

def _get_yellow_colors():
    return ["#D6B85A", "#DFC98A", "#C8A951", "#E7C27D", "#BDA55D", "#E4D00A", "#FFEF00"]

def _get_green_colors():
    return ["#3CB043", "#3A5311", "#728C69", "#AEF359", "#5DBB63", "#028A0F", "#234F1E", "#568203", "#4CBB17", "#487800"]

def _get_blue_colors():
    return ["#3944BC", "#63C5DA", "#0A1172", "#281E5D", "#1338BE", "#48AAAD", "#016064", "#2832C2", "#1F456E", "#4682B4"]

def _get_purple_colors():
    return ["#A32CC4", "#7A4988", "#601A35", "#A1045A", "#663046", "#311432", "#9867C5", "#880085"]

def _get_pink_colors():
    return ["#FC94AF", "#F25278", "#FE7D6A", "#FD5DA8", "#E3256B", "#FF6EC7"]

def _get_brown_colors():
    return ["#4B371C", "#231709", "#795C34", "#CC7722", "#65350F", "#652A0E", "#8A3324"]

def _get_black_colors():
    return ["#000000"]

def _get_grey_colors():
    return ["#808080", "#373737", "#594D5B", "#3E3D53", "#9897A9", "#63645E"]


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


def compute_scale(fig, pos, text):

    ax = fig.gca()
    master_annotation = plt.annotate(text,pos)
    bbox = master_annotation.get_window_extent()
    [x0,y0],[x1,y1] = ax.transData.inverted().transform(bbox)

    box_xlength = x1 - x0
    box_ylength = y1 - y0

    box_xlength = box_xlength * (1+3/len(text))
    box_ylength = box_ylength * 1.3

    return box_xlength, box_ylength

def compute_annotation_df(pos_list, text_list):
    annotation_df = pd.DataFrame({
        'position' : pos_list,
        'annotation' : text_list,
        'grid_coords' : np.NaN
    })
    annotation_df["grid_coords"] = annotation_df["grid_coords"].astype('object')
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
    index = grid[grid["coord"].isin(coordinates_list) == True].index
    grid.loc[index, "empty"] = False
    return grid


def find_closest_available_space(coords, grid: pd.DataFrame, x_unit, y_unit) :
    available_grid = grid.copy()
    taken_spaces = grid.query("empty == False").index
    available_grid = available_grid.drop(taken_spaces, axis= 0)
    x,y = coords
    available_grid["distance"] = np.sqrt( np.power((available_grid["x"]*x_unit - x),2) + np.power((available_grid["y"]*y_unit - y),2) )
    available_grid = available_grid.sort_values('distance').reset_index(drop= False)
    
    return available_grid.at[0, "index"]


def give_available_space(annotation_index, grid, annotation_df, x_unit, y_unit) :
    coords = annotation_df.at[annotation_index, 'position']
    space_index = find_closest_available_space(coords,grid, x_unit, y_unit)
    grid.at[space_index, "empty"] = False
    annotation_df.at[annotation_index, "grid_coords"] = grid.at[space_index, "coord"]

    return grid, annotation_df

def get_space_position(grid_coords, x_unit,y_unit) :

    x_pos = x_unit * grid_coords[0]
    y_pos = y_unit * grid_coords[1]

    return x_pos,y_pos



def write_annotation(annotation_df, x_unit, y_unit, master_length) :
    annotation_obj_list = []
    for idx in annotation_df.index :
        text = annotation_df.at[idx, "annotation"]
        grid_coords = annotation_df.at[idx, "grid_coords"]
        xy = annotation_df.at[idx,'position']
        xy_text = get_space_position(x_unit=x_unit,y_unit=y_unit, grid_coords=grid_coords)
        distance = np.sqrt(np.power((xy[0] - xy_text[0]),2) + np.power((xy[1] - xy_text[1]),2))
        xy_text = (xy_text[0] + np.random.rand() * 1.5/master_length * random_direction(), xy_text[1])        
        
        if distance >= x_unit/3 : 
            arrow_patch = {"width" : 0.25, "headwidth" : 5, 'headlength' : 2}
        else : arrow_patch = None
        annotation_obj_list.append(plt.annotate(text, xy=xy, xytext = xy_text, arrowprops= arrow_patch))

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
    master_length = len(text_list[0])
    annotation_df = compute_annotation_df(pos_list[1:],text_list[1:])
    grid = compute_grid(x_unit, y_unit)
    coords = find_grid_coordinates_list(pos_list, x_unit=x_unit, y_unit=y_unit)
    grid = fill_grid(coords, grid)

    for idx in annotation_df.index :
        grid, annotation_df = give_available_space(annotation_index= idx, annotation_df= annotation_df, grid= grid, x_unit=x_unit, y_unit=y_unit)
    annotations_obj_list = write_annotation(annotation_df,x_unit,y_unit, master_length)

    return annotations_obj_list


def is_overlapping(box1, box2) :
    if box1 == box2 : return False
    ax = plt.gca()
    [xmin1,ymin1],[xmax1,ymax1] = ax.transData.inverted().transform(box1)
    [xmin2,ymin2],[xmax2,ymax2] = ax.transData.inverted().transform(box2)

    test = 0
    if xmin1 < xmin2 and xmin2 < xmax1 : test +=1
    if ymin1 < ymin2 and ymin2 < ymax1 : test +=1
    if xmin2 < xmin1 and xmin1 < xmax2 : test +=1
    if ymin2 < ymin1 and ymin1 < ymax2 : test +=1
    return test > 1


def hide_overlapping_annotations(*annotations) :
    truth_list = []
    
    for annotation in annotations :
        test = [is_overlapping(annotation.get_window_extent(), vs_an.get_window_extent()) for vs_an in annotations]
        truth_list.append(any(test))
    
    for annotation, overlaps in zip(annotations, truth_list) :
        if overlaps :
            annotation.remove()
