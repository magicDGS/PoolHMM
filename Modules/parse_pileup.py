#!/usr/bin/python
# Modified by Daniel Gomez-Sanchez: adding opening bgzipped pileups
import sys
import getopt
import re
import os.path
import gzip

def openPileup(pileup_prefix, mode='r', log=False):
    """
    Open pileup in text format or gzipped
    pileup_prefix is the prefix of .pileup.gz or .pileup
    mode is the way of opening the file
    log=True prints if compressed or not compressed are used
    return a filehandler
    exit if IOError found
    """
    #pileup file loading
    try:
        zippedFile = pileup_prefix + ".pileup.gz"
        textFile = pileup_prefix + ".pileup"
        ## try to find the gzipped one
        if os.path.isfile(zippedFile):
            fileToOpen = zippedFile
            pileup = gzip.open(fileToOpen, mode)
        else:
            ## try to open the other one
            fileToOpen = textFile
            pileup = open(fileToOpen, mode)
        if log:
            print ("Input file: %s" % fileToOpen)
        ## return the fh
        return pileup
    except IOError:
        print ("Could not open input file %s" %  fileToOpen)
        sys.exit()


class Format:
    "Create the format of the output files"
    def __init__(self):
        self.sep = ' '
    #...

    def format(self,format,p):
        ar = []        
        if format == 'info':
            # pos chr refc nucs qual totcov eucov alleles A T C G N del valid valid_alleles derived_allele allele_code
            ar = [p['ch'], p['pos'], str(p['eucov']), p['allele_code'], p['ancestral_allele'], p['derived_allele'], str(p['removed_alleles']), str(p['unfolded']) ]
        elif format == 'freq':
	    ancestral = p['ancestral_allele']
            derived = p['derived_allele']
            nucs = p['nucs']
            ar = []
            for n in nucs:
                code = 0
		if n.upper() == ancestral:
		    ar.append(str(code))
                elif n.upper() == derived:
                    code = 1
		    ar.append(str(code))
                # ...
            #...
        elif format == 'qual':
	    ancestral = p['ancestral_allele']
            derived = p['derived_allele']
            qual = p['qual']
	    nucs = p['nucs']
            ar = []
	    i=0
            for q in qual:
		n = nucs[i]
		if n.upper() == ancestral or n.upper() == derived:
                    ar.append(str(q))
		i += 1
            #...
        #...
        return self.sep.join(ar)
    #...
#...

class Pileup_parser_ref:
    "Parse the pileup entry"
    def __init__(self,qualityEncoding,minQual,minCount,minCoverage,maxCoverage):
        self.qualityEncoding = qualityEncoding
        self.minQual = minQual
        self.minCount = minCount
        self.minCoverage = minCoverage
        self.maxCoverage = maxCoverage
        qualHash = {'illumina': lambda x : ord(x) - 64, 'sanger': lambda x : ord(x) - 33}
        if self.qualityEncoding in qualHash:
            self.encode = qualHash[self.qualityEncoding]
        else:
            print ("Quality encoding %s not recognised" % (qualityEncoding))
            sys.exit(7)
        #...
    #...

    
    def get_pileup_parser(self,line):
        if line == "\n":
            print ("pileup parser empty line provided")
            sys.exit(8)
        #...
        
        #line is splited to get ch, pos, rc, cov, nucs & qual infos
	try:
	    line_items=line.split()
        except ValueError:
            print ("could not parse pileup entry %s" % (line))
            sys.exit(9)
        #...        
	
	if len(line_items) == 6:
	    (ch, pos, rc, cov, nucs, qual) = line_items
	else:
	    print ("wrong number of columns in pileup line: %s" % (line))
	    sys.exit()

        #nucs is filtered
        #step 1
        s1 = re.findall(r"[-+](\d+)", nucs)

        for n in s1:
	        l1 = ['[-+]',n,'[ACGTNacgtn]{',n, '}']
	        pattern = r''.join(l1)
	        nucs = re.sub(pattern, '', nucs)
        #...

        #step2
        nucs = re.sub('\^.', '', nucs)

        #step3
        nucs = re.sub('\$', '', nucs)
        
        #step 4
        nucs = re.sub('\.', rc.upper(), nucs)
        
        #step 5
        nucs = re.sub(',', rc.lower(), nucs)
        #...

        #exception size of sequence != size of quality
        if len(nucs) != len(qual):
            print ("Size of sequence does not equal size of quality: %s, %s" % (nucs,line))
            sys.exit(10)
        #...
        
        #filter the pileup file by quality 
        i, ac, tc, cc, gc, nco, dell, co = 0, 0, 0, 0, 0, 0, 0, 0

        nucs_filtered = []
        qual_filtered = []
        for qc in qual:
            nc = nucs[i]
            quality = self.encode(qc)
            if quality >= self.minQual:
                co += 1
                if nc == "A" or nc == "a":
                    ac += 1
                elif nc == 'T' or nc == 't':
                    tc += 1
                elif nc == "C" or nc == "c":
                    cc += 1
                elif nc == "G" or nc == "g":
                    gc += 1
                elif nc == "N" or nc == "n":
                    nco += 1
                elif nc == "*":
                    dell += 1
                else: 
                    print ("Could not parse pileup; Unknown allele : %s in %s" % (a,line))
                    sys.exit(1)
                #...
                if nc != "*" and nc != "n" and nc != "N":
                    nucs_filtered.append(nc)
                    qual_filtered.append(quality)

            #...
            i += 1
        #...

        alar = [{'a':'A', 'c':ac}, {'a':'T', 'c':tc}, {'a':'C', 'c':cc}, {'a':'G', 'c':gc}]
        eucov = ac + tc + cc + gc
        
        if len(nucs_filtered) != len(qual_filtered):
            print ("Error: Length of nucleotides does not agree with lenght of quality string!")
            sys.exit(11)
        #...
        
        if len(nucs_filtered) != eucov:
            print ("Error : Coverage does not agree with length of nucleotides : %d n: %d " % (eucov, len(nucs_filtered)))
            sys.exit(12)
        #...
        
        # pos chr refc nucs qual totcov eucov alleles A T C G N del valid valid_alleles derived_allele allele_code
        entry={
            'pos':pos,
            'ch':ch,
            'refc':rc,
            'nucs':nucs_filtered,
            'qual':qual_filtered,
            'totcov':co,
            'eucov':eucov,
            'alleles':alar,
            'del':dell,
            'N':nco
        };

        alleles = 0
        if ac >= self.minCount:
            alleles += 1
        #...

        if tc >=  self.minCount:
            alleles += 1
        #...

        if cc >=  self.minCount:
            alleles += 1
        #...

        if gc >=  self.minCount:
            alleles += 1
        #...
        
	allele_code = "na"
        if alleles == 1:
            allele_code = "M"
        elif alleles == 2:
            allele_code = "S"
        elif alleles >2:
            allele_code = "T"
	entry['allele_code'] = allele_code        
	#..

        valid = 1
        if  entry['del'] >0 or entry['eucov'] < self.minCoverage or entry['eucov'] > self.maxCoverage:
            valid = 0
	entry["valid"] = valid

	entry['ancestral_allele']="N"
	entry['derived_allele']="N"
	entry['removed_alleles']= 0
	entry['unfolded']=1
	entry['refc']=entry['refc'].upper()

	if entry['refc'] == "A":
	    entry['ancestral_allele']="A"
	    if ac == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][0]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=ac+entry['alleles'][0]['c']
	elif entry['refc'] == "T":
	    entry['ancestral_allele']="T"
	    if tc == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][1]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=tc+entry['alleles'][0]['c']
	elif entry['refc'] == "C":
	    entry['ancestral_allele']="C"
	    if cc == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][2]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=cc+entry['alleles'][0]['c']
	elif entry['refc'] == "G":
	    entry['ancestral_allele']="G"
	    if gc == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][3]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=gc+entry['alleles'][0]['c']
	else: # entry['refc'] == 'N'
	    entry['unfolded']=0
	    entry['removed_alleles']= max(alleles-2,0)
	    if alleles >= 1:
		# sort the list alar with a bubble sort
     		alar = self.sort(alar)
		# assigns alleles			
		entry['ancestral_allele']=entry['alleles'][0]['a']
		if alleles >= 2:
        	    entry['derived_allele']=entry['alleles'][1]['a']
		    entry['eucov']=entry['alleles'][0]['c']+entry['alleles'][1]['c']

        return entry
    #...

    #Bubble sort
    def sort(self,l):
        i = 0
        changes = 1
        while (changes == 1):
            changes = 0
            i = 0
            while(i < len(l)-1):
                if l[i]['c'] < l[i+1]['c']:
                    changes = 1
                    swap = l[i+1]
                    l[i+1] = l[i]
                    l[i] = swap
                #...
                i += 1
            #...
        #...
        return l
    #..
#...

class Pileup_parser_provided:
    "Parse the pileup entry"
    def __init__(self,qualityEncoding,minQual,minCount,minCoverage,maxCoverage):
        self.qualityEncoding = qualityEncoding
        self.minQual = minQual
        self.minCount = minCount
        self.minCoverage = minCoverage
        self.maxCoverage = maxCoverage
        qualHash = {'illumina': lambda x : ord(x) - 64, 'sanger': lambda x : ord(x) - 33}
        if self.qualityEncoding in qualHash:
            self.encode = qualHash[self.qualityEncoding]
        else:
            print ("Quality encoding %s not recognised" % (qualityEncoding))
            sys.exit(7)
        #...
    #...

    def get_pileup_parser(self,line):
        if line == "\n":
            print ("pileup parser empty line provided")
            sys.exit(8)
        #...
        
        #line is splited to get ch, pos, rc, cov, nucs, qual and anc infos
        try:
	    line_items=line.split()
        except ValueError:
            print ("could not parse pileup entry %s" % (line))
            sys.exit(9)
        #...        
	
	if len(line_items) == 7:
	    (ch, pos, rc, cov, nucs, qual, anc) = line_items
	    anc=anc.upper()
	    if not ( anc == 'A' or anc == 'C' or anc == 'G' or anc == 'T' or anc == 'N'):
	    	print ("invalid ancestral allele in pileup line: %s" % (line))
	    	sys.exit()
	else:
	    print ("wrong number of columns in pileup line: %s" % (line))
	    sys.exit()
	    
        #nucs is filtered
        #step 1
        s1 = re.findall(r"[-+](\d+)", nucs)

        for n in s1:
	        l1 = ['[-+]',n,'[ACGTNacgtn]{',n, '}']
	        pattern = r''.join(l1)
	        nucs = re.sub(pattern, '', nucs)
        #...

        #step2
        nucs = re.sub('\^.', '', nucs)

        #step3
        nucs = re.sub('\$', '', nucs)
        
        #step 4
        nucs = re.sub('\.', rc.upper(), nucs)
        
        #step 5
        nucs = re.sub(',', rc.lower(), nucs)
        #...

        #exception size of sequence != size of quality
        if len(nucs) != len(qual):
            print ("Size of sequence does not equal size of quality: %s, %s" % (nucs,line))
            sys.exit(10)
        #...
        
        #filter the pileup file by quality 
        i, ac, tc, cc, gc, nco, dell, co = 0, 0, 0, 0, 0, 0, 0, 0

        nucs_filtered = []
        qual_filtered = []
        for qc in qual:
            nc = nucs[i]
            quality = self.encode(qc)
            if quality >= self.minQual:
                co += 1
                if nc == "A" or nc == "a":
                    ac += 1
                elif nc == 'T' or nc == 't':
                    tc += 1
                elif nc == "C" or nc == "c":
                    cc += 1
                elif nc == "G" or nc == "g":
                    gc += 1
                elif nc == "N" or nc == "n":
                    nco += 1
                elif nc == "*":
                    dell += 1
                else: 
                    print ("Could not parse pileup; Unknown allele : %s in %s" % (a,line))
                    sys.exit(1)
                #...
                if nc != "*" and nc != "n" and nc != "N":
                    nucs_filtered.append(nc)
                    qual_filtered.append(quality)

            #...
            i += 1
        #...

        alar = [{'a':'A', 'c':ac}, {'a':'T', 'c':tc}, {'a':'C', 'c':cc}, {'a':'G', 'c':gc}]
        eucov = ac + tc + cc + gc
        
        if len(nucs_filtered) != len(qual_filtered):
            print ("Error: Length of nucleotides does not agree with lenght of quality string!")
            sys.exit(11)
        #...
        
        if len(nucs_filtered) != eucov:
            print ("Error : Coverage does not agree with length of nucleotides : %d n: %d " % (eucov, len(nucs_filtered)))
            sys.exit(12)
        #...
        
        # pos chr refc nucs qual totcov eucov alleles A T C G N del valid valid_alleles derived_allele allele_code
        entry={
            'pos':pos,
            'ch':ch,
            'refc':rc,
            'nucs':nucs_filtered,
            'qual':qual_filtered,
            'totcov':co,
            'eucov':eucov,
            'alleles':alar,
            'del':dell,
            'N':nco
        };

        alleles = 0
        if ac >= self.minCount:
            alleles += 1
        #...

        if tc >=  self.minCount:
            alleles += 1
        #...

        if cc >=  self.minCount:
            alleles += 1
        #...

        if gc >=  self.minCount:
            alleles += 1
        #...
        
	allele_code = "na"
        if alleles == 1:
            allele_code = "M"
        elif alleles == 2:
            allele_code = "S"
        elif alleles >2:
            allele_code = "T"
	entry['allele_code'] = allele_code        
	#..

        valid = 1
        if  entry['del'] >0 or entry['eucov'] < self.minCoverage or entry['eucov'] > self.maxCoverage:
            valid = 0
	entry["valid"] = valid

	entry['ancestral_allele']="N"
	entry['derived_allele']="N"
	entry['removed_alleles']= 0
	entry['unfolded']=1
	entry['refc']=anc # here we use the provided allele instead of the reference one

	if entry['refc'] == "A":
	    entry['ancestral_allele']="A"
	    if ac == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][0]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=ac+entry['alleles'][0]['c']
	elif entry['refc'] == "T":
	    entry['ancestral_allele']="T"
	    if tc == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][1]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=tc+entry['alleles'][0]['c']
	elif entry['refc'] == "C":
	    entry['ancestral_allele']="C"
	    if cc == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][2]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=cc+entry['alleles'][0]['c']
	elif entry['refc'] == "G":
	    entry['ancestral_allele']="G"
	    if gc == 0:
		entry['removed_alleles']= max(alleles-1,0)
	    else:
		entry['removed_alleles']= max(alleles-2,0)
	    del entry['alleles'][3]
	    # sort the list alar (remaining alleles) with a bubble sort
     	    alar = self.sort(alar)
	    # assigns derived allele
	    if entry['alleles'][0]['c'] > 0:
		entry['derived_allele']=entry['alleles'][0]['a']
	    entry['eucov']=gc+entry['alleles'][0]['c']
	else: # entry['refc'] == 'N'
	    entry['unfolded']=0
	    entry['removed_alleles']= max(alleles-2,0)
	    if alleles >= 1:
		# sort the list alar with a bubble sort
     		alar = self.sort(alar)
		# assigns alleles			
		entry['ancestral_allele']=entry['alleles'][0]['a']
		if alleles >= 2:
        	    entry['derived_allele']=entry['alleles'][1]['a']
		    entry['eucov']=entry['alleles'][0]['c']+entry['alleles'][1]['c']

        return entry
    #...

    #Bubble sort
    def sort(self,l):
        i = 0
        changes = 1
        while (changes == 1):
            changes = 0
            i = 0
            while(i < len(l)-1):
                if l[i]['c'] < l[i+1]['c']:
                    changes = 1
                    swap = l[i+1]
                    l[i+1] = l[i]
                    l[i] = swap
                #...
                i += 1
            #...
        #...
        return l
    #..
#...

class Pileup_parser_folded:
    "Parse the pileup entry"
    def __init__(self,qualityEncoding,minQual,minCount,minCoverage,maxCoverage):
        self.qualityEncoding = qualityEncoding
        self.minQual = minQual
        self.minCount = minCount
        self.minCoverage = minCoverage
        self.maxCoverage = maxCoverage
        qualHash = {'illumina': lambda x : ord(x) - 64, 'sanger': lambda x : ord(x) - 33}
        if self.qualityEncoding in qualHash:
            self.encode = qualHash[self.qualityEncoding]
        else:
            print ("Quality encoding %s not recognised" % (qualityEncoding))
            sys.exit(7)
        #...
    #...

    def get_pileup_parser(self,line):
        if line == "\n":
            print ("pileup parser empty line provided")
            sys.exit(8)
        #...
        
        #line is splited to get ch, pos, rc, cov, nucs, qual and anc infos
        try:
	    line_items=line.split()
        except ValueError:
            print ("could not parse pileup entry %s" % (line))
            sys.exit(9)
        #...        
	
	if len(line_items) == 6:
	    (ch, pos, rc, cov, nucs, qual) = line_items
	else:
	    print ("wrong number of columns in pileup line: %s" % (line))
	    sys.exit()
	    
        #nucs is filtered
        #step 1
        s1 = re.findall(r"[-+](\d+)", nucs)

        for n in s1:
	        l1 = ['[-+]',n,'[ACGTNacgtn]{',n, '}']
	        pattern = r''.join(l1)
	        nucs = re.sub(pattern, '', nucs)
        #...

        #step2
        nucs = re.sub('\^.', '', nucs)

        #step3
        nucs = re.sub('\$', '', nucs)
        
        #step 4
        nucs = re.sub('\.', rc.upper(), nucs)
        
        #step 5
        nucs = re.sub(',', rc.lower(), nucs)
        #...

        #exception size of sequence != size of quality
        if len(nucs) != len(qual):
            print ("Size of sequence does not equal size of quality: %s, %s" % (nucs,line))
            sys.exit(10)
        #...
        
        #filter the pileup file by quality 
        i, ac, tc, cc, gc, nco, dell, co = 0, 0, 0, 0, 0, 0, 0, 0

        nucs_filtered = []
        qual_filtered = []
        for qc in qual:
            nc = nucs[i]
            quality = self.encode(qc)
            if quality >= self.minQual:
                co += 1
                if nc == "A" or nc == "a":
                    ac += 1
                elif nc == 'T' or nc == 't':
                    tc += 1
                elif nc == "C" or nc == "c":
                    cc += 1
                elif nc == "G" or nc == "g":
                    gc += 1
                elif nc == "N" or nc == "n":
                    nco += 1
                elif nc == "*":
                    dell += 1
                else: 
                    print ("Could not parse pileup; Unknown allele : %s in %s" % (a,line))
                    sys.exit(1)
                #...
                if nc != "*" and nc != "n" and nc != "N":
                    nucs_filtered.append(nc)
                    qual_filtered.append(quality)

            #...
            i += 1
        #...

        alar = [{'a':'A', 'c':ac}, {'a':'T', 'c':tc}, {'a':'C', 'c':cc}, {'a':'G', 'c':gc}]
        eucov = ac + tc + cc + gc
        
        if len(nucs_filtered) != len(qual_filtered):
            print ("Error: Length of nucleotides does not agree with lenght of quality string!")
            sys.exit(11)
        #...
        
        if len(nucs_filtered) != eucov:
            print ("Error : Coverage does not agree with length of nucleotides : %d n: %d " % (eucov, len(nucs_filtered)))
            sys.exit(12)
        #...
        
        # pos chr refc nucs qual totcov eucov alleles A T C G N del valid valid_alleles derived_allele allele_code
        entry={
            'pos':pos,
            'ch':ch,
            'refc':rc,
            'nucs':nucs_filtered,
            'qual':qual_filtered,
            'totcov':co,
            'eucov':eucov,
            'alleles':alar,
            'del':dell,
            'N':nco
        };

        alleles = 0
        if ac >= self.minCount:
            alleles += 1
        #...

        if tc >=  self.minCount:
            alleles += 1
        #...

        if cc >=  self.minCount:
            alleles += 1
        #...

        if gc >=  self.minCount:
            alleles += 1
        #...
        
	allele_code = "na"
        if alleles == 1:
            allele_code = "M"
        elif alleles == 2:
            allele_code = "S"
        elif alleles >2:
            allele_code = "T"
	entry['allele_code'] = allele_code        
	#..

        valid = 1
        if  entry['del'] >0 or entry['eucov'] < self.minCoverage or entry['eucov'] > self.maxCoverage:
            valid = 0
	entry["valid"] = valid

	entry['ancestral_allele']="N"
	entry['derived_allele']="N"
	entry['removed_alleles']= 0
	entry['unfolded']=0

        entry['removed_alleles']= max(alleles-2,0)
	if alleles >= 1:
	    # sort the list alar with a bubble sort
     	    alar = self.sort(alar)
	    # assigns alleles			
	    entry['ancestral_allele']=entry['alleles'][0]['a']
	    if alleles >= 2:
        	entry['derived_allele']=entry['alleles'][1]['a']
		entry['eucov']=entry['alleles'][0]['c']+entry['alleles'][1]['c']

        return entry
    #...

    #Bubble sort
    def sort(self,l):
        i = 0
        changes = 1
        while (changes == 1):
            changes = 0
            i = 0
            while(i < len(l)-1):
                if l[i]['c'] < l[i+1]['c']:
                    changes = 1
                    swap = l[i+1]
                    l[i+1] = l[i]
                    l[i] = swap
                #...
                i += 1
            #...
        #...
        return l
    #..
#...

