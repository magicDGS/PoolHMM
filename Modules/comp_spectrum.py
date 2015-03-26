import numpy as np
from EM import EM

def comp_spectrum(p,n,theta,ancestral):
    #loading infiles
    nb_sites=p.shape[0]
    # estimation of the neutral frequency spectrum
    p_neutral=EM(n,nb_sites,p,1,50,nb_sites/1000000.,theta,ancestral)
    return p_neutral
#...
