import scipy.stats as stats
import numpy as np


def ANOVA(sample_list : list, return_p_value= True, return_f_stats= False) :
    """
    Performs ANOVA test (Analysis of Variance) using `scipy.stats.f_oneway` function.
    
    Returns
    -------
    p_value

    f_stats
    """
    if len(sample_list) > 1 :
        f_stats, p_value = stats.f_oneway(*sample_list)
    else : 
        f_stats, p_value = np.NaN

    if return_p_value and return_f_stats :
        return p_value, f_stats
    elif return_p_value :
        return p_value
    elif return_f_stats :
        return f_stats