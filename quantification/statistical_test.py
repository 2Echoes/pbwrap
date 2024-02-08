import scipy.stats as stats
import numpy as np


def ANOVA(sample_list : list, return_p_value= True, return_f_stats= False, axis=0) :
    """
    Performs ANOVA test (Analysis of Variance) using `scipy.stats.f_oneway` function.
    
    Returns
    -------
    p_value

    f_stats
    """
    if len(sample_list) > 1 :
        f_stats, p_value = stats.f_oneway(*sample_list, axis=axis)
    else : 
        f_stats, p_value = np.NaN, np.NaN

    if return_p_value and return_f_stats :
        return f_stats, p_value
    elif return_p_value :
        return p_value
    elif return_f_stats :
        return f_stats


def Tukey_hsd(sample_list, return_p_value= True, return_f_stats= False):
    """
    Performs Tukey_hsd (Tukey-Kramer Honnestly significance difference) test using `scipy.stats.tukey_hsd`function.

    Returns
    -------

    """

    if len(sample_list) > 1 :
        res = stats.tukey_hsd(*sample_list)
        f_stats, p_value = res.statistic, res.pvalue
    else : 
        f_stats, p_value = np.NaN, np.NaN

    if return_p_value and return_f_stats :
        return f_stats,p_value
    elif return_p_value :
        return p_value
    elif return_f_stats :
        return f_stats

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