import numpy as np

def log_bino(n,p):
    res = np.sum(np.log(np.arange(n-p+1,n+1))) - np.sum(np.log(np.arange(1,p+1)))
    return res
#...
