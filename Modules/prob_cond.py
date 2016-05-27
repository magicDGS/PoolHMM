# Modified by Daniel Gomez-Sanchez: adding portableQueue dependency for MacOS compatibility, and opening bgzipped pileups
from portableQueue import Queue
from multiprocessing import Process, Lock
import parse_pileup as pp
import numpy as np
import time
import sys
from scipy.stats import bernoulli
from prob_cond_true_freq import prob_cond_true_freq
from comp_spectrum import comp_spectrum

#function that is parallelized
def process_probCond(qinput,qoutput,lock,pileup_prefix,parser_parameters,ratio,n,ancestral):
    pileup = pp.openPileup(pileup_prefix, 'r')
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
    
    for item in iter(qinput.get,'STOP'):
        l = []
        lock.acquire()
        pileup.seek(item[0])
        for i in range(item[1]):
            l.append(pileup.readline())
        #...
        lock.release()
        p_list = []
        for l_item in l: 
            parsed = parser.get_pileup_parser(l_item)
            if parsed['valid'] == 1:
                if bernoulli.rvs(1./ratio) == 1:
                    info = f.format('info',parsed)
                    unfolded = int(info.split()[7])
                    SE = np.fromstring(f.format('qual',parsed), dtype=float, sep=' ')
		    #print SE
                    votemp = np.fromstring(f.format('freq',parsed), dtype=int, sep=' ')
                    SEtemp = 10**(-SE/10)
                    p = prob_cond_true_freq(n,votemp,SEtemp,unfolded);
		    if np.sum(p) > 0:
                    	p_list.append(p)
                #....
            #...
        #...
        if len(p_list) != 0: #in case that all the lines parsed are not valid 
            qoutput.put(p_list)
            
        #...
    #...
    pileup.close()
#...

def prob_cond(parser_parameters,region,theta,nProcess,ratio,n,prefix,pileup_prefix,ancestral):
    lock = Lock()
    task_queue = Queue()
    done_queue = Queue()
    block = 10000
    pileup = pp.openPileup(pileup_prefix, 'rb')
   
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
    
    #for each offset except the last one
    for offset in offset_table[:-1]:
        task_queue.put([offset,block])
    #...
    
    #management of the last line_block 
    if nbLine % block != 0:
        task_queue.put([offset_table[-1],nbLine % block])
    #...
    
    del offset_table
    
    for i in range(nProcess):
        task_queue.put('STOP')
    #...
    
    for i in range(nProcess):
        p = Process(target=process_probCond,args=(task_queue,done_queue,lock,pileup_prefix,parser_parameters,ratio,n,ancestral)).start()  
    #...
    
    while task_queue.qsize() != 0:
        pass
    #...
    
    p_neutral = []
    for i in range(done_queue.qsize()):
        p_neutral += done_queue.get()
    #...
    
    p_neutral = np.array(p_neutral)
    p_neutral = comp_spectrum(p_neutral,n,theta,ancestral)
    np.savetxt(prefix + '.spectrum', np.array([p_neutral]),delimiter=' ', fmt='%.6e')
    return p_neutral
#...
