import numpy as np

def comp_emit_BF_3_SEG(n,p,p0,f_neutral,f_sel1,f_sel2,unfolded):

    # n: number of haplotypes in the pool
    # p: conditional prob of vo given z,r,SE,unfolded. size nb_sites*(n+1)
    # p0: conditional prob of 0 given z,r,SE,unfolded. size nb_sites*(n+1)
    # f_neutral: neutral frequency spectrum, size n+1
    # coeff1,coeff2: levels of hirch-hiking in the non-neutral hidden states
    # E: emission matrix, size nb_states*nb_sites
    E = np.zeros(3)
    if (unfolded==1):
        E[0] = np.sum(f_neutral * p) / (1-np.sum(f_neutral * p0))
        E[1] = np.sum(f_sel1 * p) / (1-np.sum(f_sel1 * p0))
        E[2] = np.sum(f_sel2 * p) / (1-np.sum(f_sel2 * p0))
    else:
        E[0] = np.sum(f_neutral * p) / (1-2*np.sum(f_neutral * p0))
        E[1] = np.sum(f_sel1 * p) / (1-2*np.sum(f_sel1 * p0))
        E[2] = np.sum(f_sel2 * p) / (1-2*np.sum(f_sel2 * p0))
    return E
#...
