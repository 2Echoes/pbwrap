import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter


def from_label_get_centeroidscoords(label):
    """Returns dict{"label", "centroid"}"""

    check_parameter(label = (np.ndarray))
    dim = label.ndim
    properties_dic = regionprops_table(label, properties= ["label","centroid"])
    centroid = properties_dic
    return centroid


def gene_bar_plot(rna_list: 'list[str]', values: 'list[float]', errors: 'list[float]') :
    
    is_listoflist = False
    if type(values[0]) == list : 
        if type(errors[0]) != list : raise TypeError("When passing several bar sets to plot, it is expected that several errors sets be passed as well")
        is_listoflist = True
    if type(errors[0]) == list : 
        if type(values[0]) != list : raise TypeError("When passing several errors sets to plot, it is expected that several bar sets be passed as well")
        is_listoflist = True

    color_list = ['red','blue','green','orange','purple','brown','cyan'] * (round(len(rna_list)/7) + 1)
    fig = plt.figure(figsize= (20,10))
    ax = fig.gca()
    width = 1

    if is_listoflist :
        bar_number = len(values)
        length = width/bar_number
        abs = np.arange(0,len(values[0]))
        barshift = np.arange(-(width/2 - length/2),(width/2), step = length)
        color_list = color_list[:len(barshift)]
        assert len(barshift) == len(values), "barshift : {0} ; values : {1}".format(len(barshift), len(values))
        for bar_set, error_set, shift, color in zip(values, errors, barshift, color_list) :
            X = abs - shift
            ax.bar(X, bar_set, yerr= error_set, capsize= 3, color= color, width= length, align= 'center')

    else :
        plt.bar(rna_list, values, yerr= errors, capsize= 3, color= color_list[:len(rna_list)], width= width, align= 'center')
    
    plt.axis(ymin= 0)
    
    #
    plt.xticks(range(len(rna_list)))
    xticks = ax.get_xticks()
    ax.set_xticks(xticks, labels= rna_list, rotation= 90)
    ax.axis(xmin=-0.5, xmax= len(rna_list)+0.5, ymin= 0)
    fig.subplots_adjust(bottom= 2/fig.get_size_inches()[1])

    return fig



def threshold_scatter_plot(gene_list: 'list[str]', values: 'list[float]', thresholds: 'list[float]', gene_per_col = 10) :
    #TODO
    if len(gene_list) != len(values) : raise ValueError("length of gene_list, values and thresholds must match.")
    
    fig = plt.figure(figsize=(15,15))
    gene_number = len(gene_list)
    colors = np.arange(1,1000,1000/gene_number)

    plt.scatter(thresholds, values, c= colors)
    ax = fig.gca()
    ax.axis(xmin= 0, ymin= 0)

    #Legend
    col_num = (gene_number // gene_per_col) + 1
    print("gene number ", gene_number)
    print("col num ", col_num)
    legend = ax.legend(gene_list, loc='upper left', ncols= col_num)




    return fig



def save_plot(path_output, ext):
    """Save the plot.

    Parameters
    ----------
    path_output : str
        Path to save the image (without extension).
    ext : str or List[str]
        Extension used to save the plot. If it is a list of strings, the plot
        will be saved several times.

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