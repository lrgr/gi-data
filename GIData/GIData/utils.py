import numpy as np

def sparsity(arr):
    _a = arr.reshape(-1)
    return np.isnan(_a).sum() / len(_a)