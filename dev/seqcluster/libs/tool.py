from collections import OrderedDict
import operator
import os
import copy
from libs.sam2bed import *
from libs import table,barchart,seqviz,expchart
import time

class cluster_info_obj:
    def __init__(self,clus_obj,clus_id,loci_obj,seq_obj):
        self.clus=clus_obj
        self.clusid=clus_id
        self.loci=loci_obj
        self.seq=seq_obj

class pairs:
        def __init__(self):
                self.p2={}
                self.st={}
        def set_dist(self,p,d,st):
                self.p2[p]=d
                self.st[p]=st

class sequence:
    def __init__(self, seq,freq,seq_id):
        self.seq=seq
        self.freq=freq
        self.len=len(seq)
        self.pos=[]
        self.id=seq_id
    def addpos(self,pos_id):
        self.pos.append(pos_id)

class position:
    def __init__(self,idl, chr,start,end,strand):
        self.idl=idl
        self.chr=chr
        self.start=int(start)
        self.end=int(end)
        self.strand=strand

class annotation:
    def __init__(self,db,name,strand,to5,to3):
        self.db=db
        self.name=name
        self.strand=strand
        self.to5=to5
        self.to3=to3

class dbannotation:
    def __init__(self,na):
        self.ann={}
    def adddbann(self,idl,ndba):
        self.ann[idl]=ndba
      
class cluster:
    def __init__(self, id):
        self.id=id
        self.loci={}
        self.locimax=0
        self.locilen={}
        self.dbann={}
        self.loci2seq={}
        self.idmembers={}
        self.pairs={}
        self.dist_pair={}
        self.ref=0
        self.score=0
        self.showseq=""
        self.showseq_plain=""
    def set_ref(self,r):
        self.ref=r

    def addpair(self,p1,p2,d,s):
                if not (self.pairs.has_key(p1)):
                        self.pairs[p1]=pairs()
                        self.dist_pair[p1]=1
                self.pairs[p1].set_dist(p2,d,s)
                self.dist_pair[p1]+=1
   
    def addidmember(self,ids,idl):
        lenid=len(ids)
        self.locilen[idl]=lenid
        for s in ids:
            #self.idmembers[int(c)]=0
            if not self.loci2seq.has_key(idl):
                self.loci2seq[idl]=[]
            self.loci2seq[idl].append(s)
        if (lenid>self.locimax):
            self.locimax=lenid

    def addloci(self,idl,np):
        self.loci[idl]=1
       
   
    def addmember(self,ns):
        self.members.append(ns)
       

    def adddb(self,db,ndb):
        self.dbann[db]=ndb

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

class bedaligned:
        def __init__(self,l):
            l=l.strip()
            cols = l.split("\t")
            self.chr=cols[0]
            self.start=cols[1]
            self.end=cols[2]
            self.name=cols[3]
            self.att=cols[4]
            self.strand=cols[5]

class mergealigned:
        def __init__(self,l):
                self.chr=l[0]
                self.strand=l.strand
                self.start=l.start
                self.end=l.end
                self.names=list(l.name.split(","))
                self.loci=l.score.split(",")

# def mergeclus(cur_clus_obj,clus_obj):
#     for idc in cur_clus_obj.keys():
#         #print "cluster %s " % idc
#                 #retrieve old_clus_id at clus_obj               
#         for m in cur_clus_obj[idc].idmembers.keys():
#             #print "member %s" % m          
#             for id_old in clus_obj[m].members:
#                 cur_clus_obj[idc].addmember(id_old)
   
#     return cur_clus_obj

def get_distance(p1,p2):
        if p1.strand=="+":
                d = p1.start-p2.start
        else:
                d = (p1.end-(p1.end-p1.start)) - (p2.end-(p2.end-p2.start))
        return d

def get_ini_str(n):
        #get space values
        return "".join(" " for i in range(n))

def readbowtie():
    ##read results from bowtie...no longer needed
    dict={}
    f = open("temp.map", 'r')
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        numism=line.count(">")
        mism=""
        if numism >=1:
            mism=cols[7]
       
        result="%s:%s %s %s;" % (cols[2],cols[3],cols[1],mism)
        if dict.has_key(cols[0]):
            store=dict[cols[0]]+result
            dict[cols[0]] = store  
        else:
            dict[cols[0]] = result  
    f.close()
    return dict

def what_is(file):
    f = open(file,'r')
    cols= f.readline().split("\t")
    f.close()
    if (cols[5]=="+" or cols[5]=="-"):
        return "bed"
    else:
        return "sam"

def init_numlocidb(beds):
    dict={}
    for filebed in beds:
        db= os.path.basename(filebed)
        dict[db]=0;
    return dict

def calc_complexity(nl):
    ns=len(nl)
    total=0.0
    for l in nl:
        total+=1.0/l
        #print total
    total/=ns
    return (total)

def show_seq(clus_obj,index):
    ##create vizualization of sequence along precursor for the results
    current=clus_obj.clus
    clus_seqt=clus_obj.seq
    clus_locit=clus_obj.loci
    
    itern=0
    for idc in current.keys():
        itern+=1
        #print idc
        #timestamp=str(time.time())
        timestamp=str(idc)
        seqListTemp=()
        f=open("/tmp/"+timestamp+".fa","w")
        for idl in current[idc].loci2seq.keys():
            seqListTemp=list(set(seqListTemp).union(current[idc].loci2seq[idl]))
        for s in seqListTemp:
            seq=clus_seqt[s]
            f.write(">"+s+"\n"+seq.seq+"\n")
        f.close()

        locilen_sorted=sorted(current[idc].locilen.iteritems(), key=operator.itemgetter(1),reverse=True)
        lmax=clus_locit[locilen_sorted[0][0]]
        f=open("/tmp/"+timestamp+".bed","w")
        f.write("%s\t%s\t%s\t.\t.\t%s\n" % (lmax.chr,lmax.start,lmax.end,lmax.strand))
        f.close()
        os.system("bedtools  getfasta -s -fi "+index+" -bed /tmp/"+timestamp+".bed -fo /tmp/"+timestamp+".pre.fa")
        os.system("bowtie2-build /tmp/"+timestamp+".pre.fa  /tmp/"+timestamp+".pre.ind >/dev/null 2>&1")
        os.system("bowtie2 --rdg 7,3 --mp 4 --end-to-end --no-head --no-sq  -D 20 -R 3 -N 0 -i S,1,0.8 -L 3 -f /tmp/"+timestamp+".pre.ind /tmp/"+timestamp+".fa -S /tmp/"+timestamp+".map  >>bowtie.log 2>&1")
        f=open("/tmp/"+timestamp+".map","r")
        seqpos={}
        minv=10000000
        for line in f:
            line=line.strip()
            cols=line.split("\t")
            seqpos[cols[0]]=int(cols[3])
            if minv>int(cols[3]):
                minv=int(cols[3])
        f.close()
        seqpos_sorted=sorted(seqpos.iteritems(), key=operator.itemgetter(1),reverse=False)
        showseq=""
        showseq_plain=""
        for (s,pos) in seqpos_sorted:
            ##calculate the mean expression of the sequence and change size letter
            showseq_plain+="<br><a href=javascript:loadSeq(\"%s\")>%s%s</a>" % (s,"".join("." for i in range(pos-1)),clus_seqt[s].seq)
            #showseq+=seqviz.addseq(pos-1,clus_seqt[s].len,clus_seqt[s].seq)
        #current[idc].showseq=showseq
        current[idc].showseq_plain=showseq_plain
        os.system("rm /tmp/"+timestamp+"*")
        #if (itern==1):
        #    break
    clus_obj.clus=current
    return clus_obj

def anncluster(c,clus_obj,db,type_ann):
    ##intersect positions of clusters with annotation files
    if type_ann=="gtf":
        id_sa=1
        id_ea=2
        id_sb=9
        id_eb=10
        id_id=3
        id_idl=4
        id_sta=5
        id_stb=12
        id_tag=8
    else:
        id_sa=1
        id_ea=2
        id_sb=7
        id_eb=8
        id_id=3
        id_idl=4
        id_sta=5
        id_stb=11
        id_tag=9

    clus_id=clus_obj.clus
    for cols in c.features():
        id=int(cols[id_id])
        idl=int(cols[id_idl])
        if (clus_id.has_key(id)):
            clus=clus_id[id]
            strd="-"
            sa=int(cols[id_sa])
            ea=int(cols[id_ea])
            sb=int(cols[id_sb])
            eb=int(cols[id_eb])
            if (cols[id_sta] in cols[id_stb]):
                strd="+"
            if (cols[id_sta] in "+" and cols[id_stb] in "+"):
                lento5=sa-sb+1
                lento3=ea-eb+1
            if (cols[id_sta] in "+" and cols[id_stb] in "-"):
                lento5=ea-sb+1
                lento3=sa-eb+1
            if (cols[id_sta] in "-" and cols[id_stb] in "+"):
                lento5=sa-eb+1
                lento3=ea-sb+1
            if (cols[id_sta] in "-" and cols[id_stb] in "-"):
                lento3=sa-sb+1
                lento5=ea-eb+1
            else:
                lento5=sa-sb+1
                lento3=ea-eb+1
            #print "EXISTS clus %s with db DBA %s" % (cols[9],db)
            if (clus.dbann.has_key(db)):
                #print "NEW clus %s with db DBA %s" % (cols[9],db)
                ann=annotation(db,cols[id_tag],strd,lento5,lento3)
                tdb=clus.dbann[db]
                tdb.adddbann(idl,ann)
                clus.adddb(db,tdb)
            else:
                #print "UPDATE clus %s with db DBA %s" % (cols[9],db)
                ann=annotation(db,cols[id_tag],strd,lento5,lento3)
                tdb=dbannotation(1)
                tdb.adddbann(idl,ann)
                clus.adddb(db,tdb)
           
            clus_id[id]=clus
    clus_obj.clus=clus_id
    return clus_obj

def parse_align_file(file,format):
    #parse sam files with aligned sequences
    f = open(file, 'r')
    loc_id=1
    clus_obj={}
    bedfile_clusters=""
    for line in f:
        loc_id+=1
        if format=="sam":
            line=processSAM(line)
        if line:
            a=bedaligned(line)
            #print "%s\t%s\t%s\t%s\t%s\t%s\n" % (a.chr,a.start,a.end,a.name,loc_id,a.strand)
            bedfile_clusters+="%s\t%s\t%s\t%s\t%s\t%s\n" % (a.chr,a.start,a.end,a.name,loc_id,a.strand)
       
    f.close()
    return bedfile_clusters

def parse_merge_file(c,seq_l_in,MIN_SEQ):
    #merge loci with shared sequences in same clusters
    clus_id={}
    loci_id={}
    currentClus={}
    index=0
    lindex=0
    currentClustersBed=""
    eindex=0

    for line in c.features():
        a=mergealigned(line)
        #print line
        #only keep locus with 10 or more secuences
        if (len(a.names)>=MIN_SEQ):
            exists=0
            lindex=1+lindex
            #print line
            for e in a.names:
                if (clus_id.has_key(e)):
                    
                    if exists==1 and clus_id[e]!=eindex:
                        toremove=clus_id[e]
                        #print "####already another cluster detected %s eindex %s" % (clus_id[e],eindex)
                        for e1 in currentClus[clus_id[e]].idmembers.keys():
                            clus_id[e1]=eindex
                            #move idmembers to eindex cluster
                            currentClus[eindex].idmembers[e1]=0
                            #move loci dic to eindex cluster
                            for loci_old in currentClus[toremove].loci2seq.keys():
                                if not currentClus[eindex].loci2seq.has_key(loci_old):
                                    #print "######adding lociold %s of %s to %s" % (loci_old,toremove,eindex)
                                    currentClus[eindex].addidmember(list(currentClus[toremove].loci2seq[loci_old]),loci_old)
                        #print "#####remove %s but not eindex %s" % (clus_id[e],eindex)
                        del currentClus[toremove]
                    elif exists==0:
                        exists=1
                        eindex=clus_id[e]
                    #print "exists %s in eindex %s" % (e,eindex)
                    #break;
            if exists ==1:
                #print "sequences added as %s" % eindex
                for e in a.names:
                    clus_id[e]=eindex
            else:

                index=index+1
                eindex=index
                #print "new cluster %s" % eindex
                for e in a.names:
                    clus_id[e]=eindex

            if ( not currentClus.has_key(eindex)):  
                currentClus[eindex]=cluster(eindex)

            #dev::distance from merge-start to seqs-start
            #print "%s %s" % (a.names,eindex)
            #print "%s %s %s %s %s %s" % (lindex,a.loci[0],a.names[0],a.chr,a.start,a.end)   
            newpos=position(lindex,a.chr,a.start,a.end,a.strand)
            loci_id[lindex]=newpos
            for  s in a.names:
                seq_l_in[s].addpos(lindex)
                currentClus[eindex].idmembers[s]=0
                #print s
            #currentClus[eindex].addloci(lindex)
            currentClus[eindex].addidmember(a.names,lindex)

            #dev
            #currentClustersBed=currentClustersBed + "%s\t%s\t%s\t%s\t%s\t%s\n " % (a.chr,a.start,a.end,eindex,lindex,a.strand)
    return cluster_info_obj(currentClus,clus_id,loci_id,seq_l_in)

def reduceloci(clus_obj,min_seq):
    ##reduce number of loci a cluster has
    filtered={}
    idcNew=0
    current=clus_obj.clus
    clus_idt=clus_obj.clusid
    clus_locit=clus_obj.loci
    clus_seqt=clus_obj.seq
    saved=list()
    for idc in current.keys():
        clus1=copy.deepcopy(current[idc])
        #print "idc %s " % idc
        nElements=len(clus1.loci2seq)
        #print "loci2seq %s" % clus1.loci2seq.keys()
        currentElements=nElements+1
        removeSeqs=list()
        seqListTemp=list()
        cicle=0
        seqfound=0
        while (nElements<currentElements and nElements!=0):
            cicle+=1
            #print "nC %s nE %s" % (currentElements,nElements)
            ##get loci with more number of sequences
            locilen_sorted=sorted(clus1.locilen.iteritems(), key=operator.itemgetter(1),reverse=True)
            #print "locilen %s" % locilen_sorted
            first_run=0
            ##get gigest loci
            maxseq=locilen_sorted[0][1]*1.0
            #print "max %s" % maxseq
            if maxseq>min_seq:
                for (idl,lenl) in locilen_sorted:
                    #print "starting %s %s" % (idl,lenl)
                    #print "pos %s" % clus_locit[idl].chr
                    
                    #print(clus1.loci2seq[idl])
                    #tempList=list(set(sorted(clus1.loci2seq[idl])).difference(sorted(removeSeqs)))
                    tempList=clus1.loci2seq[idl]
                    tempList=list(set(sorted(tempList)).difference(sorted(removeSeqs))) 
                    if (first_run==0):
                        seqListTemp=tempList
                    intersect=list(set(seqListTemp).intersection(tempList))
                    common=0
                    # if ("seq_389564" in tempList):
                    #     print "GOT IT:seq in idc %s cicle %s idl %s " % (idc,cicle,idl)
                    #print intersect
                    if (intersect):
                        common=len(intersect)*1.0/max(len(seqListTemp),len(tempList))
                        #print "%s %s %s %s" % (len(seqListTemp),len(tempList),
                        #len(intersect),common)
                    ##check if number of common sequences is 70% or greater than maxseq
                    if ((first_run==0 and lenl*1.0>0) or (lenl*1.0>=0.6*maxseq and common*1.0>=0.6)):
                        if (first_run==0):
                            idcNew+=1 ##creating new cluster id
                        if (not filtered.has_key(idcNew)):
                            filtered[idcNew]=cluster(idcNew)
                        #print "adding to %s" % (idcNew)
                        ##remove sequences from cluster
                        removeSeqs=list(set(tempList).union(removeSeqs))
                        # if ("seq_389564" in removeSeqs):
                        #     print "INSIDE REMOVE: seq in idc %s cicle %s idl %s newidc %s" % (idc,cicle,idl,idcNew)
                        ##adding sequences to new cluster
                        filtered[idcNew].addidmember(list(tempList),idl)
                        # if ("seq_389564" in tempList):
                        #     print "ADDED: seq in idc %s cicle %s idl %s newidc %s" % (idc,cicle,idl,idcNew)
                            # seqfound=1
                        #print  "Lt %s L2S %s IN %s C %s" % (lenl,len(clus1.loci2seq[idl]),
                        #len(intersect),common)
                        ##remove loccus from cluster
                        clus1.loci2seq.pop(idl,"None")
                        clus1.locilen.pop(idl,"None")
                    else:

                        # if ("seq_389564" in tempList):
                        #     print "###REMOVE seq in idc %s cicle %s idl %s " % (idc,cicle,idl)
                        # print "remove all sequnces in removeSeqs form current cluster"
                        ##remove all sequences in the other cluster
                        newList=list(set(sorted(tempList)).difference(sorted(removeSeqs)))                   
                        # if seqfound==1:
                        #     print "####REMOVING SEQ FROM OTHERS: %s" % idl
                        #     print ("seq_389564" in newList)
                        ##update sequences in this locus
                        clus1.locilen[idl]=len(newList)
                        clus1.loci2seq[idl]=newList
                        #print "new list %s len %s leninclus %s" % (newList,len(newList),len(clus1.loci2seq[idl]))
                        if (clus1.locilen[idl]==0):
                            clus1.loci2seq.pop(idl,"None")
                            clus1.locilen.pop(idl,"None")
                            
                    first_run=1
            else:
                ##if number of sequences < 10, just remove it
                for (idl,lenl) in locilen_sorted:
                    clus1.loci2seq.pop(idl,"None")
                    clus1.locilen.pop(idl,"None")

            nElements=len(clus1.loci2seq)
            #print(clus1.locilen.keys())
    clus_obj.clus=filtered
    return clus_obj

def generate_position_bed(clus_obj):
    ##generate file with positions in bed format
    bedaligned=""
    clus_id=clus_obj.clus
    for idc in clus_id.keys():
        clus=clus_id[idc]
        for idl in clus.loci2seq.keys():
            pos=clus_obj.loci[idl]
            bedaligned+= "%s\t%s\t%s\t%s\t%s\t%s\n" % (pos.chr,pos.start,pos.end,idc,idl,pos.strand)

    return bedaligned

def parse_ma_file(file):
    f = open(file, 'r')
    name=""
    seq_l={}
    index=1
    line=f.readline()
    line=line.strip()
    cols=line.split("\t")
    samples=cols[2:]
    for line in f:
        line=line.strip()
        cols=line.split("\t")
        exp={}
        for i in range(len(samples)):
            exp[samples[i]]=cols[i+2] 
        name=cols[0].replace(">","")
        index=index+1
        seq=cols[1]      
        new_s=sequence(seq,exp,index)
        seq_l[name]=new_s      
    f.close()
    return seq_l
