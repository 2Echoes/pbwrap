import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pbwrap.plot.bar as bar
import pbwrap.plot.utils as plot
import pbwrap.quantification.statistical_test as s_test



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

def ANOVA_test_plot(ax, data: pd.Series, significance= 0.05, colors=None, title = "ANOVA test", xlabel=None) :
    
    if len(data.index.levels) != 2 : raise ValueError("Expected 2 levels in data index : [group, treatment]. Given dim : {0}".format(len(data.index.ndim)))
    group_key = data.index.names[0]
    measure_name = data.name

    anova_p_value = data.reset_index().groupby(group_key)[measure_name].apply(s_test.ANOVA).rename('anova_p_value')

    if type(colors) != None :
        colors = pd.Series(data = list(colors), index = data.index)


    ax = bar.bar_plot(
        ax=ax,
        data= anova_p_value,
        multi_bar_plot= False,
        colors=None,
        labels= anova_p_value.index,
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

def p_value_plot(ax, data: pd.Series, statistical_test, significance= 0.05, colors=None, title = None, xlabel=None, ignore_level_test = False) :
        
    if len(data.index.levels) != 2  and not ignore_level_test: raise ValueError("Expected 2 levels in data index : [group, treatment]. Given dim : {0}".format(len(data.index.levels)))
    group_key = data.index.names[0]
    measure_name = data.name

    p_value = data.reset_index().groupby(group_key)[measure_name].apply(statistical_test).rename('p_value')
    print(p_value)

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
