import numpy as np

def hmm_viterbi(T,E,mu0):
    # pos: vector of genomic positions (in bp), nb_sites*1   
    # T: per site transition matrix of the hidden Markov process, size nb_state*nb_state
    # E: emission matrix, size nb_states*nb_sites
    # mu0: initial distribution (one should take the stationary one)
    # pred: vector of hidden states, nb_sites*1

    nb_sites = len(E)
    nb_states = np.size(T,1)
    Teff = T
    E = E.T # transpose
    Z = np.zeros((nb_states,nb_sites))
    Zarg = np.zeros((nb_states,nb_sites))

    for j in range(nb_states):
        Z[j,0]=np.log(E[j+1,0])+np.log(mu0[j])
    #...
    for i in range(1,nb_sites):
        d = int(E[0,i] - E[0,i-1])
        if d > 0:
            Teff = T**d
            Teff = np.array(Teff)
        else:
            Teff = np.identity(nb_states) * ( 1 - np.exp(-50) ) + (np.ones(nb_states) - np.identity(nb_states) ) * np.exp(-50) /2
        #...
        for j in range(nb_states):
            aux = Z[:,i-1] + np.log(Teff[:,j])
            v = np.argmax(aux) + 1
            u =np.max(aux)
            Z[j,i]= u + np.log(E[j+1,i])
            Zarg[j,i] = v
            #...
    #...
    pred = np.zeros(nb_sites)
    v = np.argmax(Z[:, nb_sites-1])+1
    pred[nb_sites-1]=v
    for i in range(nb_sites-2,-1,-1):

        pred[i] = Zarg[int(pred[i+1]-1), i+1]
    #...
    '''for i in range(len(pred)):
        pred[i] += 1
    #...'''
    return pred
#...

