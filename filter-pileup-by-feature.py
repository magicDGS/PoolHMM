#!/usr/bin/python
import sys
import collections
from optparse import OptionParser, OptionGroup


class GTFEntry:
	def __init__(self,chrom,source,feature,start,end,score,strand,offset,comment):
		self.chrom=chrom
		self.source=source
		self.feature=feature
		self.start=int(start)
		self.end=int(end)
		self.score=score
		self.strand=strand
		self.offset=offset
		self.comment=comment

class GTFReader:
	def __init__(self,file):
		self.__filename=file
		self.__filehandle=open(file,"r")

	def __iter__(self):
		return self
	
		
	@classmethod
	def readfeature(cls,file,feature):
		r=GTFReader(file)
		print "Loading gtf. file "+file
		print "Requested feature "+feature
		toret=[]
		for e in r:
			if e.feature == feature:
				toret.append(e)
		print "Finished loading; Reed " + str(len(toret))+ " features"
		return toret

	
	def next(self):
		line=""
		while(1):
			line=self.__filehandle.readline()
			if line=="":
				raise StopIteration
			line=line.rstrip('\n')
			if line != "" and line[0] != "#":
				break
		#2L	FlyBase	chromosome_band	-204333	22221	.	+	.	ID=band-21A_chromosome_band;Name=band-21A
		#2L	FlyBase	chromosome_band	-204333	-153714	.	+	.	ID=band-21A1_chromosome_band;Name=band-21A1
		e=line.split()
		e=e[0:9]
		return GTFEntry(*e)
		
		
class BinaryAnnBuilder:

	def __init__(self,annotation):
		self.__annlist=annotation
	
	def built_feature(self,chrom):
		features=self.__annlist
		
		# filter for chromosome first
		chrfeat=[]
		for feat in features:
			if feat.chrom==chrom:
				chrfeat.append(feat)
		
		# than build a binary index
		binrep=set([])
		for feature in chrfeat:	
			for i in range(feature.start,feature.end+1):
				binrep.add(i)
		return binrep


parser = OptionParser()
parser.add_option("--pileup", dest="pileup", help="A pileup file that should be filtered for a given feature")
parser.add_option("--gtf",dest="gtf",help="The annotation of the investigated species in the gtf format")
parser.add_option("--feature",dest="feature",help="The feature that should be filtered (e.g: exon, intron). Note: this feature has to be present in the gtf file")
parser.add_option("--output", dest="output", help="A filtered pileup file" )
(options, args) = parser.parse_args()

gtfentries = GTFReader.readfeature(options.gtf, options.feature)
binrep=BinaryAnnBuilder(gtfentries)

ofh=open(options.output,"w")
activeChr=""
activeRep=None

print "Starting to parse pileup file"
for line in open(options.pileup):
	line=line.rstrip()
	a=line.split("\t")
	chr=a[0]
	pos=int(a[1])
	if(chr != activeChr):
		activeChr=chr
		print "Processing chromosome "+ activeChr
		activeRep=binrep.built_feature(activeChr)
	if(pos in activeRep):
		ofh.write(line+"\n")
print "Finished"
		
	
	
