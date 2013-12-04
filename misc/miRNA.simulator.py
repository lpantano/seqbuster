from optparse import OptionParser
import sys
import os
import re
import random
from libs.objectsSim import *

def createMirna(string):
	cols = string.split("[")
	idx = 0
	obj=mirna()
	for c in cols[1:]:
		idx += 1
		c=c.replace("]","")
		p=c.split(":")
		chrom=p[0]
		p=p[1].split("-")
		start=int(p[0])
		end=int(p[1])
		
		if idx == 1:
			obj.addp5(chrom,start,end)
		else:
			obj.addp3(chrom,start,end)

	return obj
		

usagetxt = "usage: %prog  -f precurso.fa -m miRNA.str -n 10 -s hsa"
parser = OptionParser(usage=usagetxt, version="%prog 1.0")
parser.add_option("-f", "--fa", dest="fasta",
                  help="", metavar="FILE")
parser.add_option("-m", "--ma", dest="strfile",
				help="", metavar="FILE")
parser.add_option("-n", "--num", dest="numsim",
				help="", metavar="FILE")
parser.add_option("-s", "--species", dest="species",
				help="", metavar="FILE")

if len(sys.argv)<4:
    parser.print_help()
    sys.exit(1)

(options, args) = parser.parse_args()
species=options.species

pos=open(options.strfile,'r')
listmirna={}
for line in pos:	
	if line.find("[")>0:
		line = line.strip()
		name = line.split(" ")
		name = name[0].replace(">","")
		listmirna[name]=createMirna(line)
		#print name
pos.close()

print "second stage"

fas=open(options.fasta,'r')
name=""
seq=""
nt=['A','T','G','C']
for line in fas:	
	line = line.strip()
	if (line.find(">")>=0): 
		if (len(name)!=0):
			#print name
			#print len(seq)
			mir=listmirna[name].p5
			if mir!=0 and name.find(species)>=0:
				#print mir.s
				#print mir.e
				for rand in range(int(options.numsim)):
					randS=random.randint(mir.s-3,mir.s+3)
					randE=random.randint(mir.e-3,mir.e+3)
					if randS<0:
						randS=0
					if randE>mir.e:
						randE=mir.e-1					
					#print "%s %s" % (randS, randE)
					randSeq=seq[randS:randE]
					randName=name+"_"+str(randS)+":"+str(randE)
					#mutation
					isMut=random.randint(0,3)
					#print randSeq
					#print len(randSeq)
					if isMut==3:
						ntMut=random.randint(0,3)
						posMut=random.randint(0,len(randSeq)-1)
						#print posMut
						tempList=list(randSeq)
						#print tempList
						tempList[posMut]=nt[ntMut]
						randSeq=''.join(tempList)
						randName+="_mut_"+str(posMut)+":"+nt[ntMut]

					#addition
					isAdd=random.randint(0,2)
					if isAdd==2:						
						posAdd=random.randint(1,3)
						randAdd=""
						randName+="_add_"
						for numadd in range(posAdd):
							ntAdd=random.randint(0,1)
							randSeq+=nt[ntAdd]
							randName+=nt[ntAdd]

					print ">"+randName
					print randSeq
		#do random sequences
		name = line.split(" ")
		name = name[0].replace(">","")
		seq=""
	else:
		seq+=line.replace("U","T")

