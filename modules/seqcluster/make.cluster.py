import sys
import os
import re
import operator
import pybedtools
from optparse import OptionParser

from libs.tool import *

    
usagetxt = "usage: %prog  -a aligned.sam -m seq.fa -b  BED_FILE1,BED_FILE2... -o res -i index/hg19.fa"
parser = OptionParser(usage=usagetxt, version="%prog 1.0")
parser.add_option("-a", "--afile", dest="afile",
                  help="aligned file in bed/sam format", metavar="FILE")
parser.add_option("-m", "--ma", dest="ffile",
                  help="matrix file with sequences and counts for each sample", metavar="FILE")
parser.add_option("-g", "--gtf",
                   dest="gtf", help="annotate with gtf_file. It will use the 3rd column as the tag to annotate" +
                   "\nchr1    source  intergenic      1       11503   .       +       .       ",metavar="GTF")
parser.add_option("-b", "--bed",
                   dest="bed", help="annotate with bed_file. It will use the 4rd column as the tag to annotate" +
                   "\nchr1    157783  157886  snRNA   0       -",metavar="BED")
parser.add_option("-o", "--out",
                   dest="out", help="output dir",metavar="FILE")
parser.add_option("-i", "--index",
                   dest="index", help="reference fasta",metavar="FILE")

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
(options, args) = parser.parse_args()

      
if not os.path.isdir(options.out):
    print bcolors.FAIL + "the output folder doens't exists" + bcolors.ENDC
    sys.exit(1)
if not options.afile:
    print bcolors.FAIL + "-a is a mandatory parameter: alignment file in bed format" + bcolors.ENDC
    parser.print_help()
    sys.exit(1)
if not options.ffile:
    print bcolors.FAIL + "-m is a mandatory parameter: seq.ma" + bcolors.ENDC
    parser.print_help()
    sys.exit(1)
if options.bed and options.gtf:
    print bcolors.FAIL + "cannot provide -b and -g at the same time" +  bcolors.ENDC
    parser.print_help()
    sys.exit(1)


####################define variables####################    
dir_out=options.out
samplename="pronoid"

if options.bed:
    list_files=options.bed
    type_ann="bed"
if options.gtf:
    list_files=options.gtf
    type_ann="gtf"

print"output dir %s" % dir_out 
MIN_SEQ=10
db4js={}
############################################################
#try to open all files to avoid future I/O errors
try:
    f=open(options.ffile,'r')
    f.close()
    f=open(options.afile,'r')
    f.close()
    beds=list_files.split(",")
    for filebed in beds:
        f=open(filebed,'r')
        f.close()
except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit(1)
    

############read input files##################################
print "parsing matrix file"
seq_l=parse_ma_file(options.ffile)
###############################################################

#####################read aligned sequences#####################
format = what_is(options.afile)
if not format:
    print "format of aligned reads not in sam or bed"
    sys.exit(1)

print "parsing aligned file"
bed_obj=parse_align_file(options.afile,format)

###############################################################

#merge positions to create clusters: everything connected by 
#positions on the genome. 
#assumption: minimun number of common loci

a = pybedtools.BedTool(bed_obj,from_string=True)
c = a.merge(nms=True,d=20,s=True,scores="collapse")

print "Creating clusters"
clus_obj=parse_merge_file(c,seq_l,MIN_SEQ)


#####################reduce loci when possible#############################
print "Reducing multi-mapping events in the network of clusters"
setclus=reduceloci(clus_obj,MIN_SEQ,dir_out)
#sys.exit(1)
###########################################################################


#####################create sequences overview ############################
setclus=show_seq(setclus,options.index)
#sys.exit(1)
###########################################################################

#####################overlap with features#################################
print "Creating bed file"
bedfile=generate_position_bed(setclus)
a = pybedtools.BedTool(bedfile,from_string=True)
beds=[]

print "Annotating clusters"
if list_files:
    beds=list_files.split(",")
    for filebed in beds:
        db= os.path.basename(filebed)
        db4js[db]=[0,0,0]
        b = pybedtools.BedTool(filebed)
        c=a.intersect(b, wo=True)
        setclus=anncluster(c,setclus,db,type_ann) 

#############################################################

filtered=setclus.clus
clus_locit=setclus.loci
f=open(options.ffile) 
line=f.readline()
line=line.strip()
cols=line.split("\t")
samples_list=cols[2:]
f.close

#####################creating html####################
cmd="mkdir %s" % (dir_out+"/html")
print bcolors.OKBLUE+cmd+bcolors.ENDC 
os.system(cmd)
cmd="mkdir %s" % (dir_out+"/html/clusters")
print bcolors.OKBLUE+cmd+bcolors.ENDC 
os.system(cmd)
pathscript=os.path.abspath(__file__).replace("make.cluster.py","")
cmdcp="rsync -a -u --exclude='- *.' "+pathscript+"js "+pathscript+"css "+pathscript+"images "+dir_out+"/html/."
print bcolors.OKBLUE+cmdcp+bcolors.ENDC 
os.system(cmdcp)

chtml = open(dir_out+"/html/clus.html", 'w')
dbhtml = open(dir_out+"/html/annotation.html", 'w')
ccont=table.make_header("".join(map(table.make_cell_header,["id","DB"]+samples_list)))
icont=""

print "Creating outputs"
#####################creating plain text and html files#####################
out = open(dir_out+"/clus.parse.txt", 'w')
outann = open(dir_out+"/ann.tab", 'w')
outbed=dir_out+"/clus.bed"

db4js["None"]=[0,0,0]
for id in filtered.keys():
    dbsummary=""
    #if (int(id)!=0):
    clus=filtered[id]
    # "id %s" % (id)
    if (len(clus.loci2seq)>0):
        icluster_html = open(dir_out+"/html/clusters/"+str(id)+".html", 'w')
        cluster_id="C:%s" % id
        outann.write("%s\t" % cluster_id)
        ccell=table.make_cell_link(cluster_id,"clusters/%s.html" % id)
        contentDivC=table.make_div(table.make_a("<b>"+cluster_id+"</b>",cluster_id),"clus","css_id")
        out.write("C %s\n" % id)
        numloci=0
        numlocidb=init_numlocidb(beds)
        
        contentDivL=table.make_header("".join(map(table.make_cell_header,table.HEADER_L)))
        seqListTemp=()

        for idl in clus.loci2seq.keys():
            # "idl %s" % (idl)
            seqListTemp=list(set(seqListTemp).union(filtered[id].loci2seq[idl]))
            tpos=clus_locit[idl]
            numloci+=1
            out.write("L %s %s %s %s %s\n" % (tpos.idl,tpos.chr,tpos.start,tpos.end,tpos.strand))
            #outtemp.write("%s\t%s\t%s\t%s\t%s\t%s\n" %  (tpos.chr,tpos.start,tpos.end,id,0,tpos.strand))
            idl=int(tpos.idl)
            contentA=""
            hits={}
            for db in clus.dbann.keys():
                # "DBA %s" % (db)
                tdb=clus.dbann[db]
                if not hits.has_key(db):
                    hits[db]=""
                if (tdb.ann.has_key(idl)):
                    ann=tdb.ann[idl]
                    numlocidb[db]+=1
                    contentA+="%s(%s %s)-(%s,%s);" % (ann.db,ann.name,ann.strand,ann.to5,ann.to3);
                    out.write("A %s %s %s %s %s\n" % (ann.db,ann.name,ann.strand,ann.to5,ann.to3))
                    hits[db]+=ann.name+","
                    # "%s %s => %s" % (db,ann.name,hits[db])
                    #print result for consistent DB
            pos="%s:%s-%s %s" % (tpos.chr,tpos.start,tpos.end,tpos.strand)
            contentDivL+=table.make_line("".join(map(table.make_cell,[pos,contentA])))

        contentDivC+=table.make_div(table.make_table(contentDivL,"loci"),"loci","css_loci")
        #contentDivC+=table.make_div(table.make_hs_link("plainseqs"),"linkshowhide","css_seqs_link")
        contentDivC+=table.make_div(expchart.getExpDiv(),"expseqs","css_exp")
        contentDivC+=table.make_div("<pre>"+clus.showseq_plain+"</pre>","plainseqs","css_seqs")
        
        #contentDivC+=seqviz.CANVAS
        ##print seq_header
        exp=0
        seq_header="".join(map(table.make_cell_header,table.HEADER_SEQ+samples_list)) 
        contentDivS=table.make_header(seq_header)
        freq={}
        showseq=""
        allData="var allData = [{"
        for s in seqListTemp:
            allData+='"%s":[' % s
            seq=setclus.seq[s]
            out.write("S %s %s\n" % (seq.seq,seq.len))
            #showseq+="S %s %s " % (seq.seq,seq.len)
            colseqs=[seq.seq,seq.len]
            for sample in samples_list:
                allData+='{"sample":"%s","reads":"%s","color":"#FF6600"},' % (sample,int(setclus.seq[s].freq[sample]))
                colseqs.append(int(setclus.seq[s].freq[sample]))
                exp+=int(setclus.seq[s].freq[sample])
                if (not freq.has_key(sample)):
                    freq[sample]=0
                freq[sample]+=int(setclus.seq[s].freq[sample])
 
            allData=allData[:-1]+"],"
            contentDivS+=table.make_line("".join(map(table.make_cell,colseqs)))
        
        allData=allData[:-1]+"}];"

        contentDivC+=table.make_div("<pre>"+table.make_table(contentDivS,"seqs")+"</pre>","seqs","css_table")
        cons=0
        ann=0
        tmp_db4js={}
        for filebed in beds:
            db= os.path.basename(filebed)
            ratio=1.0*numlocidb[db]/numloci

            out.write("DBstat %s %s %s %s \n" % (db,numloci,numlocidb[db],ratio))
            
            if (numlocidb[db]>0):
                
                ann=1
                tmp_db4js[db]=1
                if(ratio>=0.33):
                    out.write("DB %s %s %s consistent\n" % (db,numloci,numlocidb[db]))
                    dbsummary+="DB(%s) %s/%s consistent;" % (db,numlocidb[db],numloci)
                    cons+=1
                    outann.write("%s:%s;" % (db,hits[db]))
                    
                elif (ratio<0.33):
                   
                    out.write("DB %s %s %s no-consistent\n" % (db,numloci,numlocidb[db]))
                    dbsummary+="DB(%s) %s/%s consistent;" % (db,numlocidb[db],numloci)

        for db in tmp_db4js:
            if cons==0: 
                db4js[db][2]+=exp
            if cons>1: 
                db4js[db][1]+=exp
            if cons==1: 
                db4js[db][0]+=exp

        if (ann==0): 
            db4js["None"][0]+=exp


        ccell+=table.make_cell(dbsummary)
        colseqs=[]
        outann.write("\t")
        for sample in samples_list:
            colseqs.append(freq[sample])
            outann.write("%s\t" % freq[sample])

        outann.write("\n")
        ccell+="".join(map(table.make_cell,colseqs))
        ccont+=table.make_line(ccell)
        icont=table.make_div(contentDivC,cluster_id,"css_clus")
        #idiv=make_div(clus,clus,"cluster")
    icont=table.make_html(icont,expchart.addgraph(allData),samplename)         
    icluster_html.write(icont)
    icluster_html.close()

listbar=[]
for db in db4js.keys():
    listbar.append({"args":db,"unique":db4js[db][0],"multiple":db4js[db][1],"unconsistently":db4js[db][2]})
listkeys=["unique","multiple","unconsistently"]

dbhtml.write(barchart.createhtml(listbar,listkeys))
dbhtml.close()
ccont=table.make_table(ccont,samplename)
ccont=table.make_jshtml(ccont,samplename)         
chtml.write(ccont)
chtml.close()
out.close()
outann.close()

############################################################

####################data for genome browser####################
#cmd2intersect="mergeBed -s -nms -i %s > %s" % (temporalfile,temporalbedfile)
#print bcolors.OKBLUE+cmd2intersect+bcolors.ENDC
#os.system(cmd2intersect)

#header="track name=%(name)s description='%(name)s' " % {'name':samplename}
#cmd2intersect="echo '%s';awk '{print $1,$2,$3,$4,$5,$5 }'  %s | sed 's/ /\t/g'> %s " % (header,temporalbedfile,outbed)
#print bcolors.OKBLUE+cmd2intersect+bcolors.ENDC
#os.system(cmd2intersect)

############################################################





