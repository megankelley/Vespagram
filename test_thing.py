from itertools import izip_longest
import numpy as np

def sum_n_lists(LoL, stride):
    args = [iter(LoL)] * stride
    
    sums = []
    for it in izip_longest(*args, fillvalue=None):
        np_arrays = [np.asarray(l) for l in it]
        
        np_sum = np.zeros(len(np_arrays[0]))
        for a in np_arrays:
            np_sum += a
        
        sums.append(np_sum)
    
    return sums


a = [[1, 2, 3],
     [4, 5, 6],
     [7, 8, 9],
     [1, 2, 3]]

print sum_n_lists(a, 2)
