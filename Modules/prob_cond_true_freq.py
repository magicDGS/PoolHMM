import numpy as np
from prob_cond_true_freq_unfolded import prob_cond_true_freq_unfolded

def prob_cond_true_freq(n,vo,SE,unfolded):
    #vo: vector of alleles (1 or 0), size r
    #SE: vector of sequencing error probabilities, size r
    #unfolded: 1 if unfolded, 0 else.
    #p: conditional prob of vo given z,r,SE,unfolded. size n+1*1
    o = sum(vo)
    r = len(vo)
    aux = np.ones(r,'int')
    if (unfolded==1):
        p = prob_cond_true_freq_unfolded(n,r,vo,SE)
    else:
        p = (prob_cond_true_freq_unfolded(n,r,vo,SE) + prob_cond_true_freq_unfolded(n,r,aux-vo,SE))/2
    #...
    return p
#...

