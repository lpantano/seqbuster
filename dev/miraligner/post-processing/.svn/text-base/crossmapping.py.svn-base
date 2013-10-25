#!/usr/bin/python
import sys
import os
from optparse import OptionParser
from time import clock, time


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def readbowtie(namefile):
    dict={}
    f = open(namefile, 'r')
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

        

usagetxt = "usage: %prog  -f MIRALIGNER_OUTPUT -a  BOWTIE_INDEX"
parser = OptionParser(usage=usagetxt, version="%prog 1.0")
parser.add_option("-f", "--file", dest="filename",
                  help="output from miraligner", metavar="FILE")
parser.add_option("-a", "--align",
                   dest="align", help="map sequences with bowtie using this genome index",metavar="INDEX")
parser.add_option("-b", "--bed",
                   dest="mapped", help="bed format: sequences already mapped on to a genome, available in next update",metavar="mapped")
parser.add_option('-t', action="store_true", default=False)

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
    
(options, args) = parser.parse_args()

      
if not options.filename:
    print bcolors.FAIL + "-f is a mandatory parameter: output from miraligner" + bcolors.ENDC
    parser.print_help()
    sys.exit(1)


bowtie={}
prefix=str(time())
abspath= os.path.abspath(options.filename)
path=abspath.replace(options.filename,"")
if (path==""):
    path="./"


dirout=options.filename+".results/"

if options.t:
    cmd="mkdir %s" % (dirout)
    print bcolors.OKBLUE+cmd+bcolors.ENDC 
    os.system(cmd)


pathscript=os.path.abspath(__file__).replace("crossmapping.py","")
header=""
if options.align:

    tempfa=path+prefix+".fa"
    tempmap=path+prefix+".map"
    #tempfa="kk.fa"
    #tempmap="kk.map"
     #print path
    #sys.exit(0)
    print bcolors.OKGREEN+"selected -a parameter"+bcolors.ENDC
    f = open(options.filename, 'r')
    fao = open(tempfa, 'w')
    print "%sreading file %s %s" % (bcolors.OKBLUE,options.filename,bcolors.ENDC) 
    header=f.readline()
    header=header.strip()
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        #print cols[0]
        fao.write(">%s\n%s\n" % (cols[0],cols[0]))
    fao.close()

    #execute bowtie
    cmd="bowtie -f -v 1 -a --best --strata -M 3 %s %s %s" % (options.align,tempfa,tempmap)
    print bcolors.OKBLUE+cmd+bcolors.ENDC 
    os.system(cmd)
    
    f.close()
    fao.close()
    #open bowtie output
    bowtie=readbowtie(tempmap)
    #bowtie=readbowtie("1332426993.7.map")


f = open(options.filename, 'r')
print "%sreading file %s %s" % (bcolors.OKBLUE,options.filename,bcolors.ENDC) 
header=f.readline()
out = open(options.filename+".cm", 'w')
if options.t:
    outhtml = open(dirout+"result.html", 'w')
    outputhtml.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\"><html><head><meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" /><title>miralignet output</title><style type=\"text/css\" title=\"currentStyle\">@import \"css/demo_page.css\";	@import \"css/jquery.dataTables.css\";</style><script type=\"text/javascript\" language=\"javascript\" src=\"js/jquery.js\"></script><script type=\"text/javascript\" language=\"javascript\" src=\"js/jquery.dataTables.js\"></script><script type=\"text/javascript\" charset=\"utf-8\">$(document).ready(function() {$('#mirna').dataTable();} );</script></head><body id=\"dt_example\">")
    outhtml.write("<table cellpadding=\"0\" cellspacing=\"0\" border=\"0\" class=\"display\" id=\"mirna\" width=\"300ex\">\n")
    outhtml.write("<thead><tr><th>id</th><th>freq/name</th><th>miR</th><th>start</th><th>end</th><th>substitution</th><th>5 trimming</th><th>3 trimming</th><th>3 addition</th><th>5 flanking region</th><th>3 flanking region</th><th>DB</th><th>ambiguity</th><th>hits</th><th>state</tr></thead><tbody>")

out.write(header+"\tstatus\n")
for line in f:
    line = line.strip()
    cols = line.split("\t")
    #print line
    linehtml=line.replace("\t","</td><td>")
    #linehtml='</td><td>'.join(cols)
    if bowtie.has_key(cols[0]):
        
        numism=int(bowtie[cols[0]].count(">"))
        numloc=int(bowtie[cols[0]].count(";"))

        if cols[5]!="0" and numism==0:
            out.write("%s\t%s\tproblematic\n" % (line,bowtie[cols[0]]))
            if options.t:
                outhtml.write("<tr><td>%s</td><td>%s</td><td>better match on genome w/o mismatch</td></tr>\n" % (linehtml,bowtie[cols[0]]))
        elif cols[5]!="0" and numism/numloc>=1.0:
            out.write("%s\t%s\tconfirmed\n" % (line,bowtie[cols[0]]))
            if options.t:
                outhtml.write("<tr><td>%s</td><td>%s</td><td>confirmed</td></tr>\n" % (linehtml,bowtie[cols[0]]))
        elif cols[5]!="0" and numism/numloc<1.0:
            out.write("%s\t%s\tconfirmed\n" % (line,bowtie[cols[0]]))
            if options.t:
                outhtml.write("<tr><td>%s</td><td>%s</td><td>confirmed</td></tr>\n" % (linehtml,bowtie[cols[0]]))
        elif cols[6]!=0 and numism/numloc<(len(cols[6])-1):
            out.write("%s\t%s\tproblematic\n" % (line,bowtie[cols[0]]))
            if options.t:
                outhtml.write("<tr><td>%s</td><td>%s</td><td>better match on genome w/o addition</td></tr>\n" % (linehtml,bowtie[cols[0]]))
        else:
            out.write("%s\t%s\tconfirmed\n" % (line,bowtie[cols[0]]))
            if options.t:
                outhtml.write("<tr><td>%s</td><td>%s</td><td>confirmed</td></tr>\n" % (linehtml,bowtie[cols[0]]))
    
    else:
         out.write("%s\tNA\tconfirmed\n" % line)
         if options.t:
             outhtml.write("<tr><td>%s</td><td>NA</td><td>confirmed</td></tr>\n" % (linehtml))

if options.t:
    outhtml.write("</tbody></table>")
    outhtml.close()

out.close()
f.close()  

cmdrm="rm "+path+prefix+"*"
#print bcolors.OKBLUE+cmdrm+bcolors.ENDC 
os.system(cmdrm)  

if options.t:
    cmdcp="rsync -a -u --exclude='- *.' "+pathscript+"js "+pathscript+"css "+pathscript+"images "+dirout+"."
    print bcolors.OKBLUE+cmdcp+bcolors.ENDC 
    os.system(cmdcp)

print "%s Finish cross-mapping annotation. Your outfile file is %s folder %s" % (bcolors.OKBLUE,options.filename,bcolors.ENDC) 



