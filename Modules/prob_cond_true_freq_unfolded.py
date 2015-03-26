import numpy as np

def prob_cond_true_freq_unfolded(n,r,vo,SE):
# n: number of haplotypes in the pool
# r: number of reads at the position
# vo: vector of alleles (1 if derived, 0 otherwise), size r
# SE: vector of sequencing error probabilities, size r
# p: conditional prob of vo given z,r,SE,unfolded. size n+1*1
    p = np.zeros(n+1)
    for i in range(r):
        if (vo[i] == 1):
            
            mat = np.log((((1-SE[i])*np.arange(n+1.))/n) + (SE[i]*(1-np.arange(n+1.)/n)))   
        else:
            mat = np.log(((SE[i]*np.arange(n+1.)/n))+((1-SE[i])*(1-np.arange(n+1.)/n)))
        #...
        p += mat
    #...
    p = np.exp(p)
    return p
#...
