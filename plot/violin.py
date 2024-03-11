"""
This submodules groups all function related to violin plots making from base plot to result plots.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pbwrap.utils import check_parameter
from itertools import chain
from pbwrap.plot.utils import get_colors_list


def violin_plot(
        ax: plt.Axes, distributions, labels=None, sub_labels=None, colors=None, xlabel=None, ylabel=None, title=None, y_axis= None,
        linewith = 2, line_color = 'black', alpha= 0.6,
        vertical_plot=True, showmedians= True, showextrema = True, showmeans=False, mean_size = 45, multi_violin_plot= False
                ) :

    """
    Basis for violin plots.
    """

    #Check parameters

    check_parameter(ax = plt.Axes, xlabel = (type(None), str), ylabel = (type(None), str), vertical_plot= bool, title = (type(None),str), linewith = (int, float), line_color= str, y_axis= (type(None),tuple,list), multi_violin_plot= bool)
    if multi_violin_plot :
        for distrib in distributions :
            check_parameter(distrib = (list, tuple, np.ndarray, pd.DataFrame, pd.Series))

    if type(labels) == type(None) : pass
    elif len(distributions) != len(labels) : raise ValueError("Length of labels and distributions must match.")

    if type(colors) == type(None): pass
    elif len(distributions) != len(colors) and not multi_violin_plot: raise ValueError("Length of colors and distributions must match.")
    
    if type(y_axis) == type(None) : pass
    elif len(y_axis) != 2 : raise ValueError("y_axis is expected to be of length 2 : (ymin, ymax).")
    
    if multi_violin_plot :
        if type(colors) != type(None) :
            if len(colors) == len(distributions) :
                new_colors = []
                for color,distrib in zip(colors, distributions) :
                    new_colors.extend([color]*len(distrib))
                colors = new_colors

        max_individual_violin_number = max([len(distrib) for distrib in distributions]) + 1 #Is the maximum number of violin plotted for one set.
        positions, ticks_positions = multi_violin_plot_positions(distributions)
        distributions = list(chain(*distributions))
        assert len(distributions) == len(positions), "AssertionError : multi_distributions wrongly flattened : positions : {0}, distributions {1}".format(len(positions), len(distributions))

        #colors
        if type(colors) != type(None) :
            if len(colors) == len(distributions) : pass
            else : raise ValueError("Length of colors must either match length of distributions or the number of element in distributions")
    
    else :
        positions = np.arange(1, len(distributions) + 1)
        ticks_positions = np.arange(1, len(distributions) + 1)
        max_individual_violin_number = 1
    #Plot
    violin_plot = ax.violinplot(
        distributions,
        positions=positions, 
        vert= vertical_plot, 
        showmedians=showmedians,
        showextrema= showextrema,
        )

    if type(labels) == type(None) :
        labels = np.arange(1, len(distributions) + 1)
    
    xticks = ax.set_xticks(ticks_positions, labels=labels)


    ax.set_xlim(0.25, len(labels) * max_individual_violin_number + 0.75)
    if type(colors) == type(None) : colors = get_colors_list(len(violin_plot['bodies']))
    for violin, color in zip(violin_plot['bodies'], colors) :
        violin.set_facecolor(color)
        violin.set_alpha(alpha)

    for collection_name in ['cbars', 'cmins', 'cmaxes', 'cmedians'] :
        collection = violin_plot.get(collection_name)
        if type(collection) != type(None) :    
            collection.set_color(line_color)
            collection.set_linewidth(linewith)
    
    if type(y_axis) != type(None) :
        axis = list(ax.axis())
        if y_axis[0] != None : axis[2] = y_axis[0]
        if y_axis[1] != None : axis[3] = y_axis[1]
        ax.axis(axis)
    else : axis = ax.axis()

    if type(sub_labels) != type(None) :
        if len(sub_labels) != len(distributions) : raise ValueError("Length of sub_labels must match number of violins to plot.")
        for sub_label, x_position in zip(sub_labels, positions) :
            ax.text(x = x_position, y= axis[2], s= sub_label, ha= 'center')
    
    if showmeans :
        means = [np.mean(distrib) for distrib in distributions]
        ax.scatter(positions, means, c= colors, s= mean_size, linewidths=0.5, edgecolors='black')

    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)
    if type(ylabel) != type(None) : ax.set_ylabel(ylabel)
    if type(title) != type(None) : ax.set_title(title)

    return ax


def multi_violin_plot_positions(distributions) :

    max_individual_violin_number = max([len(distrib) for distrib in distributions]) + 1#Is the maximum number of violin plotted for one set.

    positions = []
    ticks_positions = []
    for distrib_number, distrib in enumerate(distributions) :
        positions.extend(list(
            np.arange(1, len(distrib) + 1) + (distrib_number * max_individual_violin_number) if len(distrib) > 1 
            else [distrib_number * max_individual_violin_number + (max_individual_violin_number-1)/2 + 1]
        ))

        ticks_positions.append(
            distrib_number * max_individual_violin_number + (len(distrib)-1)/2 + 1 if len(distrib) > 1
            else distrib_number * max_individual_violin_number + (max_individual_violin_number-1)/2 + 1
        )

    return positions, ticks_positions


# #test
# plt.figure()
# data = [
#     [np.random.rand(100) for i in range(np.random.randint(1,5))]
#     for i in range(5)]
# for num, d in enumerate(data) :
#     print("set {0} : {1} violins".format(num, len(d)))

# sub_lab = sum([['{0} sublabel'.format(len(i))]*len(i) for i in data], [])
# print(sub_lab)

# ax = plt.gca()
# violin_plot(ax, data, labels= ['A','B','C','D','E'], colors= ['black', 'red', 'blue', 'orange', 'green'], sub_labels=sub_lab, multi_violin_plot=True,
#             xlabel= 'Label', ylabel= 'Y label', title= 'TITLE')
# plt.show()