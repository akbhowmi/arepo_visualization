import numpy as np
import sys


def apply_cuts(index,data,inverse=False):

    if inverse:
        newdata = [tup[np.in1d(range(len(tup)),index)==0] for tup in data]
    else:
        newdata = [tup[index] for tup in data]
    return newdata
