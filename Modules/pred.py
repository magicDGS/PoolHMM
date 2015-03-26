from cree_transit_3 import cree_transit_3
from scipy.linalg import expm3
from hmm_posterior import hmm_posterior
from hmm_viterbi import hmm_viterbi
from prop_pred_SEG_3 import prop_pred_SEG_3
import numpy as np

def prediction(k,prefix):
    # parameters of the HMM
    pn=0.25
    ps=0.25
    A = cree_transit_3(pn,ps,k,0)
    T = np.matrix(expm3(A))
    mu0 = np.array([pn,1-ps-pn,ps])
    E = np.loadtxt(prefix + '.segemit', delimiter = ' ', usecols = (1,2,3,4))
    print ('HMM parameters loaded')

    # computation of the posterior
    post = hmm_posterior(T,E,mu0)
    post_file = open(prefix + '.post','w')
    for i in range(len(post)):
	post_file.write( str(int(E[i,0])) + ' ' + str(post[i]) + '\n')
    post_file.close()
    print ('Posterior probabilities computed')

    # prediction of the hidden states
    pred = hmm_viterbi(T,E,mu0)
    pred_file = open(prefix + '.pred','w')
    for i in range(len(pred)):
	pred_file.write( str(int(E[i,0])) + ' ' + str(int(pred[i])) + '\n')
    pred_file.close()

    stat = prop_pred_SEG_3(E[:,0],pred,post)
    stat_file = open(prefix + '.stat','w')
    stat_file.write( str(int(stat[0][1])) + '\n')
    for i in range(1,len(stat)):
	stat_file.write( str(int(stat[i][0])) + ' ' + str(int(stat[i][1])) + ' ' + str(stat[i][2]) + '\n')
    stat_file.close()

    print ('Prediction of hidden states finished')
    if stat[0][0] == 0:
	print ('No selective sweep was found')
#...
