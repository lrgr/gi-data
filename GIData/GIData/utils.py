import numpy as np

def sparsity(arr):
    _a = arr.reshape(-1)
    return np.isnan(_a).sum() / len(_a)

def multi_index_duplicated(multi_index, key):
    ''' returns returns tuples in multi index that are duplicated, w.r.t to the given key '''
    return multi_index[multi_index.get_level_values(key).duplicated(keep=False)].values
