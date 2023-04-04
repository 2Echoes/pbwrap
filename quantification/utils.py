import numpy as np


def unzip(lis) :
    """from a list of coordinates return a list of list where each rows has all coordinate from an axis"""

    res = []
    for dim in zip(*lis):
        res += [list(dim)]
    print(res)
    return res