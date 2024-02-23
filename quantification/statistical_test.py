import scipy.stats as stats
import pingouin
import pandas as pd
import numpy as np

#Note : I've set the min_pvalue according to float decimal precesion being around 17 decimals.

def alexandergovern(sample_list : list, return_p_value= True, return_f_stats= False, axis=0, min_pvalue = 10e-19) :
    """
    Performs ANOVA-like test (for sample with unequal variance) using `scipy.stats.alexandergovern` function.
    
    Returns
    -------
    p_value

    f_stats
    """
    if len(sample_list) > 1 :
        res = stats.alexandergovern(*sample_list, nan_policy='raise')
        f_stats, p_value = res.statistic, res.pvalue
        if p_value < min_pvalue : p_value = min_pvalue
    else : 
        f_stats, p_value = np.NaN, np.NaN

    if return_p_value and return_f_stats :
        return f_stats, p_value
    elif return_p_value :
        return p_value
    elif return_f_stats :
        return f_stats

def ANOVA(sample_list : list, return_p_value= True, return_f_stats= False, axis=0, min_pvalue = 10e-19) :
    """
    Performs ANOVA test (Analysis of Variance) using `scipy.stats.f_oneway` function.
    
    Returns
    -------
    p_value

    f_stats
    """
    if len(sample_list) > 1 :
        f_stats, p_value = stats.f_oneway(*sample_list, axis=axis)
        if p_value < min_pvalue : p_value = min_pvalue
    else : 
        f_stats, p_value = np.NaN, np.NaN

    if return_p_value and return_f_stats :
        return f_stats, p_value
    elif return_p_value :
        return p_value
    elif return_f_stats :
        return f_stats


def Tukey_hsd(sample_list, return_p_value= True, return_f_stats= False, min_p_value= 10e-19):
    """
    Performs Tukey_hsd (Tukey-Kramer Honnestly significance difference) test using `scipy.stats.tukey_hsd`function.

    Returns
    -------

    """

    if len(sample_list) > 1 :
        res = stats.tukey_hsd(*sample_list)
        f_stats, p_value = res.statistic, res.pvalue
        p_value[p_value < min_p_value] = min_p_value
    else : 
        f_stats, p_value = np.NaN, np.NaN

    if return_p_value and return_f_stats :
        return f_stats,p_value
    elif return_p_value :
        return p_value
    elif return_f_stats :
        return f_stats

def games_howell(sample_list, min_p_value= 10e-19) :
    
    sample_number = len(sample_list)

    if sample_number > 1 :
        sample_ids = []
        for sample_idx, sample in enumerate(sample_list) :
            sample_ids += [sample_idx] * len(sample)

        df = pd.DataFrame({
            'index' : sample_ids,
            'measure' : np.concatenate(sample_list)
        })

        res = pingouin.pairwise_gameshowell(df, dv= 'measure', between= 'index').loc[:,['A', 'B', 'pval']]
        p_value = np.identity(sample_number)
        for sample_idx in res['A'].unique() :
            p_value[sample_idx, (sample_idx+1):] =  res[res['A'] == sample_idx].loc[:,'pval']
            p_value[(sample_idx+1):, sample_idx] =  res[res['A'] == sample_idx].loc[:,'pval']
        
        p_value[p_value < min_p_value] = min_p_value
    
    else : 
        p_value = np.NaN, np.NaN

    return p_value

def chi_squared_test(sample_list, return_p_value= True, return_statics= False) :
    """
    Performs chi-squared test that the categorical data has even frequencies (null hypothesis).
    Test performed with `scipy.stats.chisquare` 
    """
    print(sample_list)
    res = [stats.chisquare(sample,f_exp= None) for sample in sample_list]
    print('res : \n',res)
    statistic, p_value = zip(*res)

    if return_p_value and return_statics :
        return statistic, p_value
    elif return_p_value :
        return p_value
    elif return_statics :
        statistic