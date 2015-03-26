import numpy as np
from log_bino import log_bino

def proba_nielsen(n,p_neutral,coeff):
    # n: number of haplotypes
    # p_neutral: probabilities for frequencies from 0 to n
    # coeff: intensity of hichhiking
    # p_sel: probabilities for frequencies from 0 to n

    p = 1 - np.exp(-coeff)

    # number of lineages escaping the sweep

    p_e = np.zeros(n+1) 

    for k in range(n+1):
        p_e[k] = np.exp(log_bino(n,k)) * p**k * (1-p)**(n-k)
    #...  

    # probabilities of observing j mutants in a neutral sample of size h
    mat_p = np.zeros((n+1,n+1))
    for h in range(1,n+1):
        #cpt=h
        for j in range(h+1):
            for i in range(j,n-h+j+1):          
                aux = log_bino(i,j) + log_bino(n-i,h-j) - log_bino(n,h)
                mat_p[j,h] +=  p_neutral[i] * np.exp(aux)
            #...
        #...
    #...

    #probabilities after selection
    p_sel = np.zeros(n+1)
    for b in range(n+1):
        p_sel[b] = p_neutral[b] * p_e[n]
        for k in range(max(b-1,0),n):
            p_sel[b]=p_sel[b]+p_e[k]*mat_p[b,k+1]*(k+1-b)/(k+1)
        #...
        for k in range(min(n-b,n-1),n):
            p_sel[b]=p_sel[b]+p_e[k]*mat_p[b+1-n+k,k+1]*(b+1-n+k)/(k+1)
        #...
    #...
    return p_sel

