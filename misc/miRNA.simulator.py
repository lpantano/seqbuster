from optparse import OptionParser
import sys
import os
import re
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
		

usagetxt = "usage: %prog  -f precurso.fa -m miRNA.str -n 10"
parser = OptionParser(usage=usagetxt, version="%prog 1.0")
parser.add_option("-f", "--fa", dest="fasta",
                  help="", metavar="FILE")
parser.add_option("-m", "--ma", dest="strfile",
				help="", metavar="FILE")
parser.add_option("-n", "--num", dest="numsim",
				help="", metavar="FILE")

if len(sys.argv)<4:
    parser.print_help()
    sys.exit(1)

(options, args) = parser.parse_args()

pos=open(options.strfile,'r')
listmirna={}
for line in pos:	
	if line.find("[")>0:
		line = line.strip()
		name = line.split(" ")
		name = name[0].replace(">","")
		listmirna[name]=createMirna(line)
		print name
pos.close()

print "second stage"

fas=open(options.fasta,'r')
name=""
seq=""
for line in fas:	
	line = line.strip()
	if line.find(">")>=0:
		print name		
		print seq
		#do random sequences
		name = line.split(" ")
		name = name[0].replace(">","")
		seq=""
	else:
		seq+=line

