import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import regionprops_table
from bigfish.stack import check_parameter
from sklearn import linear_model



def from_label_get_centeroidscoords(label):
    """Returns dict{"label", "centroid"}"""

    check_parameter(label = (np.ndarray))
    dim = label.ndim
    properties_dic = regionprops_table(label, properties= ["label","centroid"])
    centroid = properties_dic
    return centroid


def gene_bar_plot(rna_list: 'list[str]', values: 'list[float]', errors: 'list[float]', legend: 'list[str]'= None, width = 0.8, error_width = 3) :

    #Exception thrower
    is_listoflist = False

    #multi data set
    if type(values[0]) == list : 
        if type(errors[0]) != list : raise TypeError("When passing several bar sets to plot, it is expected that several errors sets be passed as well")
        is_listoflist = True
    if type(errors[0]) == list : 
        if type(values[0]) != list : raise TypeError("When passing several errors sets to plot, it is expected that several bar sets be passed as well")
        is_listoflist = True

    #len list matches
    if not is_listoflist :
        if not(len(rna_list) == len(values) == len(errors)) : raise ValueError("rna, values and errors lengths must match")
    else :
        #Set lengths match
        if not(len(values) == len(errors)) and legend == None : raise ValueError("value sets and errors sets lengths must match")
        elif not(len(values) == len(errors) == len(legend)) : raise ValueError("value sets and errors sets and legend lengths must match")
        #Data lengths match
        for set in range(0,len(values)) :
            if not(len(rna_list) == len(values[set]) == len(errors[set])) : raise ValueError("values and errors lengths must match")

    
    #Init plot
    color_list = ['red','blue','green','orange','purple','brown','cyan'] * (round(len(rna_list)/7) + 1)
    fig = plt.figure(figsize= (20,10))
    ax = fig.gca()

    #Case when several bars are plotted for each genes
    if is_listoflist :
        bar_number = len(values)
        length = width/bar_number
        error_width /= bar_number
        abs = np.arange(0,len(values[0]))
        barshift = np.arange(-(width/2 - length/2),(width/2), step = length)
        color_list = color_list[:len(barshift)]
        assert len(barshift) == len(values), "barshift : {0} ; values : {1}".format(len(barshift), len(values))
        for bar_set, error_set, shift, color, label in zip(values, errors, barshift, color_list, legend) :
            X = abs - shift
            if legend != None : ax.bar(X, bar_set, yerr= error_set, capsize= error_width, color= color, width= length, align= 'center', label= label)
            else : ax.bar(X, bar_set, yerr= error_set, capsize= error_width, color= color, width= length, align= 'center')
        if legend != None : ax.legend()
    
    #Case one bar per gene
    else :
        plt.bar(rna_list, values, yerr= errors, capsize= error_width, color= color_list[:len(rna_list)], width= width, align= 'center')
    
    plt.axis(ymin= 0)
    
    #Spacing
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
        


def get_simple_linear_regression(X: np.array, Y: np.array) :
    X = np.array(X).reshape(-1,1)
    Y = np.array(Y)
    lin_model = linear_model.LinearRegression()
    lin_model.fit(X,Y)

    return lin_model.coef_[0], lin_model.intercept_