# Modified by Daniel Gomez-Sanchez: adding portableQueue dependency for MacOS compatibility, and opening bgzipped pileups
from portableQueue import Queue
import numpy as np
from prob_cond_true_freq import prob_cond_true_freq
from multiprocessing import Process, Lock
import parse_pileup as pp

class A:
    "class A"
    def __init__(self):
        self.list_chro = []
        self.list_pos = []
	self.list_anc = []
	self.list_der = []
        self.list_u = []
    #...
#...

#function that sorts a list of A object by the list_pos[0]
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
def process_estim(qinput,qoutput,lock,parser_parameters,pileup_prefix,n,p_neutral,ancestral):
    print 'process starts'

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
        estim_tab = A()
        for l_item in l:
            parsed =parser.get_pileup_parser(l_item)
            if parsed['valid'] == 1:
                info = f.format('info',parsed) 
                unfolded = int(info.split()[7])
                SE = np.fromstring(f.format('qual',parsed) , dtype=float, sep=' ')
                votemp = np.fromstring(f.format('freq',parsed), dtype=int, sep=' ')
                SEtemp = 10**(-SE/10)
                estim_tab.list_chro.append(info.split()[0])
                estim_tab.list_pos.append(int(info.split()[1]))
		estim_tab.list_anc.append(info.split()[4])
		estim_tab.list_der.append(info.split()[5])
		if unfolded == 1:
                    estim_tab.list_u.append(np.argmax(p_neutral * prob_cond_true_freq(n,votemp,SEtemp,1)))
		else:
		    estim_tab.list_u.append(np.argmax((p_neutral + p_neutral[::-1])* prob_cond_true_freq(n,votemp,SEtemp,1)))
		    Ltemp=len(estim_tab.list_u)-1
		    if estim_tab.list_u[Ltemp]>(n/2):
			estim_tab.list_u[Ltemp]=n-estim_tab.list_u[Ltemp]
			estim_tab.list_anc[Ltemp]=info.split()[5]
			estim_tab.list_der[Ltemp]=info.split()[4]
            #...
        #...
        #print '+1'
        if len(estim_tab.list_u) != 0: #in case that all the lines parsed are not valid
            qoutput.put(estim_tab)
        #...
    #...
    print 'process stops'
    pileup.close()
#...

def estimation(parser_parameters, region, nProcess, n, prefix,p_neutral,pileup_prefix,ancestral):
    lock = Lock()
    task_queue = Queue()
    done_queue = Queue()
    block = 10000
    pileup = pp.openPileup(pileup_prefix, 'r')

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
    
    #for each offset expept the last one
    for offset in offset_table[:-1]:
        task_queue.put([offset,block])
    #...
    
    #management of the last line_block 
    if nbLine % block != 0:
        task_queue.put([offset_table[-1],nbLine % block])
    #...
    
    for i in range(nProcess):
        task_queue.put('STOP')
    #...
    
    for i in range(nProcess):
        p = Process(target=process_estim,args=(task_queue, done_queue,lock,parser_parameters,pileup_prefix,n,p_neutral,ancestral)).start()
    #...
    
    while task_queue.qsize() != 0:
        pass
    #...
    
    estim = []
    for i in range(done_queue.qsize()):
        estim.append(done_queue.get())
    #...
    if len(estim) > 1:
        quickSort(estim)
    #...
    fic_estim = open(prefix + '.estim', 'w')
    for item in estim:
        for i in range(len(item.list_pos)):
            fic_estim.write( item.list_chro[i] + ' ' + str(item.list_pos[i]) + ' ' + item.list_anc[i] + ' ' + item.list_der[i] + ' ' + str(item.list_u[i]) + '\n')
        #...
    #...
    fic_estim.close()
#...
