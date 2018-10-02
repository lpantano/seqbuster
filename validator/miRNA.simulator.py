from __future__ import print_function
from optparse import OptionParser
import sys
import os
import re
import random
import numpy

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
                start=int(p[0]) - 1
                end=int(p[1])
                if idx == 1:
                        obj.addp5(chrom,start,end)
                else:
                        obj.addp3(chrom,start,end)
        return obj

class pos:
	def __init__(self,n,s,e):
		self.s=s
		self.e=e
		self.id=n

class mirna:
	def __init__(self):
		self.p5=0
		self.p3=0
	def addp5(self,n,s,e):
		#print("adding p5")
		self.p5=pos(n,s,e)
		#print(self.p5)
	def addp3(self,n,s,e):
		self.p3=pos(n,s,e)

usagetxt = "usage: %prog  -f precurso.fa -m miRNA.str -n 50 -s hsa -e"

parser = OptionParser(usage=usagetxt, version="%prog 1.0")
parser.add_option("-f", "--fa", dest="fasta",
                  help="", metavar="FILE")
parser.add_option("-m", "--ma", dest="strfile",
                                help="", metavar="FILE")
parser.add_option("-n", "--num", dest="numsim",
                                help="")
parser.add_option("-s", "--species", dest="species",
                                help="", metavar="FILE")
parser.add_option("-e", "--exp", dest="exp", action="store_true",
                                help="give expression", default=False)
parser.add_option("-p", "--prefix", help="output name")


if len(sys.argv)<4:
    parser.print_help()
    sys.exit(1)

(options, args) = parser.parse_args()
species = options.species
out_fa = "%s.fa" % options.prefix
out_fq = "%s.fq" % options.prefix

pos_mirna = open(options.strfile, 'r')
listmirna = {}
for line in pos_mirna:
    if line.find("[") > 0:
        line = line.strip()
        name = line.split(" ")
        name = name[0].replace(">", "")
        listmirna[name] = createMirna(line)
pos_mirna.close()

data = {}
fas = open(options.fasta, 'r')
fout = open(out_fa, 'w')
fqout = open(out_fq, 'w')
name = ""
seq = ""
nt = ['A', 'T', 'G', 'C']
for line in fas:
        line = line.strip()
        if (line.find(">")>=0):
                if (len(name)!=0):
                        mir=listmirna[name].p5
                        if mir!=0 and name.find(species)>=0:
                                for rand in range(int(options.numsim)):
                                        randS=random.randint(mir.s-2,mir.s+2)
                                        randE=random.randint(mir.e-2,mir.e+2)
                                        if randS<1:
                                                randS=1
                                        if randE>len(seq):
                                                randE=mir.e-1
                                        randSeq=seq[randS-1:randE]
                                        randName=name+"_"+mir.id+"_"+str(randS)+":"+str(randE)
                                        randName=randName+"_"+str(randS-mir.s-1)+":"+str(randE-mir.e)
                                        isMut=random.randint(0,3)
                                        mut_tag = "_mut:null"
                                        if isMut==3:
                                                ntMut=random.randint(0,3)
                                                posMut=random.randint(0,len(randSeq)-1)
                                                tempList=list(randSeq)
                                                if tempList[posMut] == nt[ntMut]:
                                                    ntMut -= 1
                                                    if ntMut < 0 :
                                                        ntMut +=2
                                                tempList[posMut]=nt[ntMut]
                                                randSeq=''.join(tempList)
                                                mut_tag="_mut:"+str(posMut+1)+nt[ntMut]
                                        randName+=mut_tag
                                        isAdd=random.randint(0,2)
                                        add_tag = "_add:null"
                                        if isAdd==2:
                                                posAdd=random.randint(1,3)
                                                randAdd=""
                                                add_tag="_add:"
                                                for numadd in range(posAdd):
                                                        ntAdd=random.randint(0,1)
                                                        randSeq+=nt[ntAdd]
                                                        add_tag+=nt[ntAdd]
                                        randName+=add_tag
                                        exp = ""
                                        if not data.has_key(randSeq):
                                            if options.exp:
                                                trial =random.randint(1, 100)
                                                p = random.randint(1, 50) / 50.0
                                                exp = "_x%s" % numpy.random.negative_binomial(trial, p, 1)[0]
                                            print("@%s%s" % (randName, exp), file=fqout, end="\n")
                                            print("%s" % randSeq, file=fqout, end="\n")
                                            print("+\n%s" % "".join(["I"] * len(randSeq)), file=fqout, end="\n")
                                            print(">%s%s" % (randName, exp), file=fout, end="\n")
                                            print("%s" % randSeq, file=fout, end="\n")
                                            data[randSeq]=1
                name = line.split(" ")
                name = name[0].replace(">", "")
                seq = ""
        else:
                seq += line.replace("U", "T")
