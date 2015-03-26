import numpy as np

def prop_pred_SEG_3(pos,pred,post):
    # returns useful statistics about the predicted sequence in segregating sites format
    nb_sweep = 0
    lim_sweep = []
    lim_sweep.append([0,0,0])
    L = len(pred)
    sweep = 0
    ind_min = 0
    ind_max = L-1

    for i in range (0,L):
        if pred[i] == 3:
            if sweep == 0:
                sweep = 1
                nb_sweep += 1
                ind_min = i
            #...
        else:
            if sweep == 1:
                ind_max =i - 1
                lim_sweep.append([pos[ind_min],pos[ind_max],max(-np.log10(1-post[ind_min:ind_max]))])
		sweep = 0
            #...    
        #...
    #...
    if sweep == 1:
        ind_max = L-1
        lim_sweep.append([pos[ind_min],pos[ind_max],max(-np.log10(1-post[ind_min:ind_max]))])
    #...
    lim_sweep[0]=[min(nb_sweep,1),nb_sweep,0]
    return lim_sweep
#...
