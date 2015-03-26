import numpy as np

def hmm_posterior(T,E,mu0):
    # pos: vector of genomic positions (in bp), nb_sites*1   
    # T: per site transition matrix of the hidden Markov process, size nb_state*nb_state
    # E: emission matrix, size nb_states*nb_sites
    # mu0: initial distribution (one should take the stationary one)
    # post: vector of posterior prob, nb_sites*1

    nb_sites = len(E)
    nb_states = np.size(T,1)
    Teff = T
    E = E.T
    alpha = np.zeros((nb_states,nb_sites))
    beta = np.ones((nb_states,nb_sites))  
    # algo forward
    for j in  range(nb_states):
        alpha[j,0] = E[j+1,0] * mu0[j]
    #...
    
    alpha[:,0] /=  np.sum(alpha[:, 0])
    for i in range(1, nb_sites):
        d= int(E[0,i]-E[0,i-1])
        if d > 0:
            Teff = T**d
            Teff = np.array(Teff)
        else:
            Teff = np.identity(nb_states) * ( 1 - np.exp(-50) ) + (np.ones(nb_states) - np.identity(nb_states) ) * np.exp(-50) /2
        #...
        for j in range(nb_states):
            alpha[j,i] = E[j+1,i] * np.sum(Teff[:, j] * alpha[:, i-1])
        #...
        alpha[:,i] /= np.sum(alpha[:, i])
    #...

    # algo backward
    for i in range(nb_sites-2,-1,-1):
        d = int(E[0,i+1] - E[0,i])
        if d > 0:
            Teff = T**d
            Teff = np.array(Teff)
        else:
            Teff = np.identity(nb_states)
            #...
        for j in range(0,nb_states):                             
            beta[j,i] = np.sum(E[1:,i+1] * Teff[j].T * beta[:,i+1])
            #...
        beta[:,i] /= np.sum(beta[:,i])
    #...
    # combination
    post = np.zeros(nb_sites)
    for i in range(nb_sites):
        v = alpha[:,i] * beta[:,i]
        v = v / np.sum(v)
        post[i] = v[2]
    #...
    return post
#...
