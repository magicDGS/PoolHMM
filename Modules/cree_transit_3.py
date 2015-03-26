import numpy as np

def cree_transit_3(ps,pn,fact,epsilon):
    # epsilon must be small compared to fact
    A = np.zeros((3,3))
    A[0,1] = fact
    A[0,2] = epsilon
    A[0,0] = -A[0,1] - A[0,2]
    q = 2/(1-ps-pn) * ((fact+epsilon) * pn-epsilon*ps)
    A[1,0] = q/2
    A[1,2] = q/2
    A[1,1] = -A[1,0] - A[1,2]
    r = (q*(1-ps-pn)-fact*pn)/ps
    A[2,0] = epsilon
    A[2,1] = r
    A[2,2] = -A[2,1] - A[2,0]
    return A
#...
