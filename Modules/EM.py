import numpy as np

def EM(n,nb_sites,p,nb_dep,maxiter,maxerr,theta,ancestral):
    # n: number of haplotypes in the pool
    # r: nb of reads, size nb_sites*1
    # p: conditional prob of vo given z,r,SE,unfolded. size nb_sites*(n+1)
    # nb_dep: number of independent starting points of the EM
    # maxiter: maximum number of iterations of the EM for one given start
    # maxerr: largest distance tolerated between two consecutive likelihoods
    # v_new: frequency spectrum, size n+1

    lik_best = -nb_sites**2
    v_best = np.zeros(n+1)
    for j in range(nb_dep):

        #simulation of the EM starting point from the expected neutral spectrum
        v_old= np.zeros(n+1)
        v_old[1:]= theta / (np.arange(1,n+1))
        v_old[0] = 1 - sum(v_old[1:])

	if ancestral == 'unknown':
	    v_old = (v_old + v_old[::-1])/2

        lik_old2= -nb_sites*10000
        iter=0
        err=1000
        #iterations
        while ((err>maxerr) and (iter < maxiter)): 
            lik_old = 0
            v_new = np.zeros(n+1)
            for i in range(nb_sites):
                aux_old = v_old * p[i]
                sumAux_old = np.sum(aux_old)
                lik_old +=  np.log(sumAux_old)
                v_new += (aux_old/sumAux_old)
              #...
            v_new /= nb_sites          
            err = lik_old - lik_old2
            iter += 1
            v_old = v_new
            lik_old2 = lik_old
        #...

        if (lik_old > lik_best):
            v_best = v_new
            lik_best = lik_old
        #...
    #...
    v_new = v_best
    return v_new
#...

