# Modified by Daniel Gomez-Sanchez: adding portableQueue dependency for MacOS compatibility, and opening bgzipped pileups
from portableQueue import Queue
from multiprocessing import Process, Lock
import parse_pileup as pp
import numpy as np
from prob_cond_true_freq import prob_cond_true_freq
from comp_emit_BF_3_SEG import comp_emit_BF_3_SEG
from proba_nielsen import proba_nielsen
import time

class A:
    "class A"
    def __init__(self):
        self.list_chro = []
        self.list_pos = []
        self.list_E = []
    #...
#...

#function that sort a list of A object by the list_pos[0]
def quickSort(L):
 
    def trirap(L, g, d):
        pivot = L[(g+d)//2].list_pos[0]
        i = g
        j = d
        while True:
            while L[i].list_pos[0]<pivot:
                i+=1
            while L[j].list_pos[0]>pivot:
                j-=1
            if i>j:
                break
            if i<j:
                L[i], L[j] = L[j], L[i]
            i+=1
            j-=1
        if g<j:
            trirap(L,g,j)
        if i<d:
            trirap(L,i,d)
 
    g=0
    d=len(L)-1
    trirap(L,g,d)
#...

#function that is parallelized
def procress_emit(qinput,qoutput,lock,pileup_prefix,parser_parameters,n,p_neutral,f_sel1,f_sel2,ancestral):
    qualityEncoding = parser_parameters[0]
    minQual = parser_parameters[1]
    minCount = parser_parameters[2]
    minCoverage = parser_parameters[3]
    maxCoverage = parser_parameters[4]
    
    #creation of the parser object
    if ancestral == "provided":
	parser = pp.Pileup_parser_provided(qualityEncoding,minQual,minCount,minCoverage,maxCoverage)
    elif ancestral == "unknown":
	parser = pp.Pileup_parser_folded(qualityEncoding,minQual,minCount,minCoverage,maxCoverage)
    else:
	parser = pp.Pileup_parser_ref(qualityEncoding,minQual,minCount,minCoverage,maxCoverage)
    f = pp.Format()
    
    pileup = pp.openPileup(pileup_prefix, 'r')
    for item in iter(qinput.get,'STOP'):
        l = []
        pileup.seek(item[0])
        for i in range(item[1]):
            l.append(pileup.readline())
        #...
        p = A()
        for l_item in l:
            parsed = parser.get_pileup_parser(l_item)      
            if parsed['valid'] == 1:
                info = f.format('info',parsed)
                if info.split()[5] != "N":                
                    unfolded = int(info.split()[7])
                    votemp = np.fromstring(f.format('freq',parsed), dtype=int, sep=' ')
                    SE = np.fromstring(f.format('qual',parsed), dtype=float, sep=' ')
                    SEtemp = 10**(-SE/10)
                    pseg = prob_cond_true_freq(n,votemp,SEtemp,unfolded)
		    votemp = np.zeros(len(votemp))
                    p0 = prob_cond_true_freq(n,votemp,SEtemp,unfolded)
		    E = comp_emit_BF_3_SEG(n,pseg,p0,p_neutral,f_sel1,f_sel2,unfolded)
		    if E[0] > 0 and E[1] > 0 and E[2] > 0:                  	
                    	p.list_chro.append(info.split()[0])
                    	p.list_pos.append(int(info.split()[1]))
                    	p.list_E.append(E)
                #...
            #...
        #...
        if len(p.list_E) != 0: #in case that all the lines parsed are not valid 
            qoutput.put(p)
        #...
    #...
    pileup.close()
#...

def comp_emit_seg_direct(parser_parameters,region,nProcess,n,prefix,p_neutral,pileup_prefix,ancestral):
    
    emit = open(prefix + '.segemit','w')
    lock = Lock()
    task_queue = Queue()
    done_queue = Queue()
    block = 10000
    pileup = pp.openPileup(pileup_prefix, 'rb')

    # computation of cond prob for segsites
    if region:
        chro = region[0]
        start = region[1]
        end = region[2]
        offset_default = pileup.tell()
        pileup_line = pileup.readline()
        a = pileup_line.split()[0]
        while(a != chro):
            offset_default = pileup.tell()
            pileup_line = pileup.readline()
            try:
                a = pileup_line.split()[0]
            except IndexError:
                #if the pileup_line can't be splited, that's the end of the file
                print ('ERROR : chro %s not found' % (chro))
                sys.exit()
            #...
        #...
        if start:
            a = int(pileup_line.split()[1])
            b = pileup_line.split()[0]
            if a > end:
                print ('ERROR : interval\'s positions not found.')
                sys.exit()
            #...
            while a < start and b == chro:
                offset_default = pileup.tell()
                pileup_line = pileup.readline()
                try:
                    a = int(pileup_line.split()[1])
                    b = pileup_line.split()[0]
                except IndexError:
                    #if the pileup_line can't be splited, that's the end of the file
                    print ('ERROR : interval\'s positions not found.')
                #...
            #...
            if b != chro:
                print b
                print ('ERROR : interval\'s positions not found.')
                sys.exit()
        #...
        offset_table = [offset_default]
        nbLine = 0
        split_pileup = pileup_line.split()
        while split_pileup[0] == chro:
            if start:
                if int(split_pileup[1]) > end:
                    break
                #...
             #...
            nbLine += 1
            if nbLine % block == 0:
                offset_table.append(pileup.tell())
            #...
            pileup_line = pileup.readline()
            split_pileup = pileup_line.split()
            if len(split_pileup) == 0:
                break
            #... 
        #...
    #...
    else:
        offset_table = [0]
        nbLine = 0
        pileup_line = pileup.readline()
        while(pileup_line != ''): #if pileup_line == '', that's the end of the file
            nbLine += 1
            if nbLine % block == 0:
                offset_table.append(pileup.tell())
            #...
            pileup_line = pileup.readline()
        #...
    #...
    pileup.close()
    coeff1 = 0.7
    coeff2 = 0.2
    f_sel1 = proba_nielsen(n,p_neutral,coeff1)
    print ('f_sel1 loaded')
    f_sel2 = proba_nielsen(n,p_neutral,coeff2)
    print ('f_sel2 loaded')
    
    for i in range(nProcess):
        p = Process(target=procress_emit,args=(task_queue, done_queue,lock,pileup_prefix,parser_parameters,n,p_neutral,f_sel1,f_sel2,ancestral)).start()
    #...
    
    #for each offset expept the last one
    for offset in offset_table[:-1]:
        task_queue.put([offset,block])
    #...
    
    #management of the last line_block 
    if nbLine % block != 0:
        task_queue.put([offset_table[-1],nbLine % block])
    #...
    del offset_table
    del f_sel1
    del f_sel2
    
    for i in range(nProcess):
        task_queue.put('STOP')
    #...
    while task_queue.qsize() != 0:
        pass
    #...
    mat = []
    for i in range(done_queue.qsize()):
        mat.append(done_queue.get())
    #...
    if len(mat) > 1:
        quickSort(mat)
    #...

    for item in mat:
        for i in range(len(item.list_pos)):
            emit.write( item.list_chro[i] + ' ' + str(item.list_pos[i]) + ' ' + str(item.list_E[i][0]) + ' ' + str(item.list_E[i][1]) + ' ' + str(item.list_E[i][2]) + '\n')
        #...
    #...
    emit.close()
#...
