import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pbwrap.plot.bar as bar
import pbwrap.plot.utils as plot
import pbwrap.quantification.statistical_test as s_test

from matplotlib.colors import SymLogNorm



def variance_plot(ax:plt.Axes, group_list: 'pd.Series[list]', labels= None, colors=None, normalise_var = True, index_normalisation=None) :

    if type(labels) == type(None) : labels = np.arange(len(group_list))

    if not type(group_list) == pd.Series :
        group_list = pd.Series(group_list, labels= labels)

    if normalise_var :
        ylabel = 'Normalised sample variance'
    else : 
        ylabel = 'Sample variance'

    #var computation

    var_data = pd.Series(dtype= float)
    if normalise_var :
        if not type(index_normalisation) in (list, np.array, tuple) : index_normalisation = np.zeros(len(group_list), dtype= int)
        if len(index_normalisation) != len(group_list) : index_normalisation = np.zeros(len(group_list))
        normal_constants = [np.var(group_list.at[sample_index][norm_index]) for sample_index, norm_index in zip(group_list.index, index_normalisation)]

        for index, norm_constant in zip(group_list.index, normal_constants) :
            var_data.at[index] = [np.var(sample) / norm_constant for sample in group_list.at[index]]
    
    else :
        for index in group_list.index :
            var_data.at[index] = [np.var(sample) for sample in group_list.at[index]]

    ax = bar.bar_plot(
        ax=ax,
        data= var_data,
        labels= labels,
        colors= colors,
        multi_bar_plot= True,
        ylabel= ylabel

    )

    return ax

def ANOVA_Tukey_hsd(ax, data: pd.Series, h0_treatment, significance= 0.05, colors=None, title = "ANOVA test followed by Tukey honnestly significance difference test", xlabel=None) :

    if len(data.index.levels) != 2 : raise ValueError("Expected 2 levels in data index : [group, treatment]. Given dim : {0}".format(len(data.index.ndim)))
    group_key, treatment_key = data.index.names[0], data.index.names[1]
    measure_name = data.name

    anova_p_value = data.reset_index().groupby(group_key)[measure_name].apply(s_test.ANOVA).rename('anova_p_value')
    rna_anova_valid = anova_p_value[anova_p_value < significance].index

    if type(colors) != None :
        colors = pd.Series(data = list(colors), index = data.index).loc[rna_anova_valid]

    treatment_list = data.loc[rna_anova_valid].reset_index().groupby(group_key)[treatment_key].apply(list) #nom des traitements après cutoff

    tukey_p_value = data.loc[rna_anova_valid].reset_index().groupby(group_key)[measure_name].apply(s_test.Tukey_hsd).rename("tukey_hsd_p_value")
    
    h0_treatment_index = pd.Series(dtype= float)
    for index in treatment_list.index :
        h0_treatment_index.at[index] = treatment_list.at[index].index(h0_treatment)

    for index in tukey_p_value.index :
        position = h0_treatment_index.at[index]
        tukey_p_value.at[index] = tukey_p_value.at[index][position]


    ax = bar.bar_plot(
        ax=ax,
        data= tukey_p_value,
        multi_bar_plot= True,
        colors=colors,
        labels= tukey_p_value.index,
        ylabel= "p-value",
        title= title,
        xlabel= xlabel,
        logscale= True
    )

    plot.plot_horizontal_bar(significance)

    return ax

def ANOVA_test_plot(ax, data: pd.Series, significance= 0.05, title = "ANOVA test", xlabel=None, xtick_label = None, group=False) :
    
    if group :
        group_key = data.index.names[0]
        measure_name = data.name
        anova_p_value = data.reset_index().groupby(group_key)[measure_name].apply(s_test.ANOVA).rename('anova_p_value')
    else :
        anova_p_value = [s_test.ANOVA(data)]
    
    if type(xtick_label) == type(None) :
        xtick_label = anova_p_value.index
        
    ax = bar.bar_plot(
        ax=ax,
        data= anova_p_value,
        multi_bar_plot= False,
        colors=None,
        labels= xtick_label,
        ylabel= "p-value",
        title= title,
        xlabel= xlabel,
        logscale= True
    )

    plot.plot_horizontal_bar(significance)

    return ax

def alexandergovern_test_plot(ax, data: pd.Series, significance= 0.05, title = "ANOVA test", xlabel=None, xtick_label = None, group=False) :
    """
    
    """
    
    if group :
        group_key = data.index.names[0]
        measure_name = data.name
        p_value = data.reset_index().groupby(group_key)[measure_name].apply(s_test.alexandergovern).rename('p_value')
    else :
        p_value = [s_test.alexandergovern(data)]
    
    if type(xtick_label) == type(None) and group :
        xtick_label = p_value.index
    elif type(xtick_label) == type(None) and not group :
        xtick_label = ['sample']
        
    ax = bar.bar_plot(
        ax=ax,
        data= p_value,
        multi_bar_plot= False,
        colors=None,
        labels= xtick_label,
        ylabel= "p-value",
        title= title,
        xlabel= xlabel,
        logscale= True
    )

    plot.plot_horizontal_bar(significance)

    return ax
    
def Tukey_hsd_plot(ax, data: pd.Series, h0_treatment, significance= 0.05, colors=None, title = "Tukey honnestly significance difference test", xlabel=None) :

    if len(data.index.levels) != 2 : raise ValueError("Expected 2 levels in data index : [group, treatment]. Given dim : {0}".format(len(data.index.ndim)))
    group_key, treatment_key = data.index.names[0], data.index.names[1]
    measure_name = data.name

    if type(colors) != None :
        colors = pd.Series(data = list(colors), index = data.index)

    treatment_list = data.reset_index().groupby(group_key)[treatment_key].apply(list) #nom des traitements après cutoff

    tukey_p_value = data.reset_index().groupby(group_key)[measure_name].apply(s_test.Tukey_hsd).rename("tukey_hsd_p_value")
    drop_index = tukey_p_value[tukey_p_value.isna()].index
    tukey_p_value = tukey_p_value.drop(drop_index, axis= 0)
    colors = colors.drop(drop_index, axis= 0)
    
    h0_treatment_index = pd.Series(dtype= float)
    for index in treatment_list.index :
        h0_treatment_index.at[index] = treatment_list.at[index].index(h0_treatment) if h0_treatment in treatment_list.at[index] else 0

    for index in tukey_p_value.index :
        position = h0_treatment_index.at[index]
        tukey_p_value.at[index] = tukey_p_value.at[index][position]


    ax = bar.bar_plot(
        ax=ax,
        data= tukey_p_value,
        multi_bar_plot= True,
        colors=colors,
        labels= tukey_p_value.index,
        ylabel= "p-value",
        title= title,
        xlabel= xlabel,
        logscale= True
    )

    plot.plot_horizontal_bar(significance)

    return ax

def pairwise_stest_tile(ax: plt.Axes, data: pd.DataFrame, measure=None, groupby_key=None, significance = 0.01, title= 'pairwise p-values',xlabel=None, test= 'tukey-hsd') :
    """
    Expected : Df with groupby_key as the treatment key to group on and measure the measure to compute pairwise pvalues on.

    Parameters
    ----------
        ax : matplotlib.pyplot.Axes
            axes to plot on
        data : pd.DataFrame/pd.Series with multi-index
        measure : str
            key to measure to compute p-values on should be passed only if data is a DataFrame that needs to be grouped.
        groupb_key : str
            key to group data on should be passed only if data is a DataFrame that needs to be grouped.
        significance : float
            Significance of the test : used for colorscale, should be between 0 and 1
        title : str
        test : str
            statistical test to perform : 'tukey-hsd' or 'gameshowell'

            > tukey-hsd : Tukey-kramer Honnestly Significant differences appropriate for samples with equal or unequal size but same finite variance. 
            Performed with `scipy.stats.tukey_hsd`
            
            > gameshowell : Pairwise test with correction for unequal sample size and unequal variance. 
            Performed with `pingouin.pairwise_gameshowell`.
    """

    if (type(measure) != type(None) and type(groupby_key) == type(None)) or (type(measure) == type(None) and type(groupby_key) != type(None)):
        raise ValueError("groupby_key and measure should be passed or set to None together.")
    elif type(measure) == type(None) and type(data) != pd.Series :
        raise ValueError("if groupby_key and measure set to None data should be a pd.Series with index corresponding to group") 

    if type(groupby_key) != type(None) and type(measure) == type(None) :
        data = data.groupby(groupby_key)[measure].apply(list)
    
    treatment_list = list(data.index)
    if test == 'tukey-hsd' : tukey_p_value = np.array(s_test.Tukey_hsd(data))
    elif test == 'gameshowell' : tukey_p_value = s_test.games_howell(data)
    else : raise ValueError('test paremeter should be either "tukey-hsd" or "gameshowell".')

    #plot
    mesh = ax.pcolormesh(
        np.flip(tukey_p_value, axis=0),
        edgecolors= 'black',
        norm= SymLogNorm(vmin= 10e-19, vmax= 1, linthresh= significance, linscale= -np.log10(significance)),
        cmap= 'bwr',
        )

    ticks_positions = np.arange(len(treatment_list)) + 0.5
    ax.set_xticks(ticks_positions, treatment_list)
    ax.set_yticks(ticks_positions, treatment_list[::-1])
    ax.xaxis.set_ticks_position('top')
    
    #color bar
    cbar = plt.colorbar(mesh, ax=ax, location= 'right', ticks = [1, significance, 1e-1, 1e-3, 1e-18])


    if type(title) != type(None) : ax.set_title(title)
    if type(xlabel) != type(None) : ax.set_xlabel(xlabel)

    return ax
    
def p_value_plot(ax : plt.Axes, data: pd.Series, statistical_test, significance= 0.05, colors=None, title = None, xlabel=None, ignore_level_test = False) :
        
    if len(data.index.levels) != 2  and not ignore_level_test: raise ValueError("Expected 2 levels in data index : [group, treatment]. Given dim : {0}".format(len(data.index.levels)))
    group_key = data.index.names[0]
    measure_name = data.name

    p_value = data.reset_index().groupby(group_key)[measure_name].apply(statistical_test).rename('p_value')

    if type(colors) != None :
        colors = pd.Series(data = list(colors), index = data.index)

    ax = bar.bar_plot(
        ax=ax,
        data= p_value,
        multi_bar_plot= False,
        colors=None,
        labels= p_value.index,
        ylabel= "p-value",
        title= title,
        xlabel= xlabel,
        logscale= True
    )

    plot.plot_horizontal_bar(significance)

    return ax
