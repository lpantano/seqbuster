library(isomiRs)


#setwd("path2seqbuster/R/isomiR_package/test")

#seq path to source code
files<-c("y0d2.hsa.fa.ad.new.mirna",
          "y0d34.hsa.fa.ad.new.mirna",
         "y0d204.hsa.fa.ad.new.mirna",
         "y66d0.hsa.fa.ad.new.mirna",
         "y80d0.hsa.fa.ad.new.mirna",
         "y88d0.hsa.fa.ad.new.mirna"
         )


d<-data.frame(condition=c("nb","nb","nb","o","o","o"))
row.names(d)<-paste(d[,1],1:3,sep="")

obj<-loadIso(files,d,skip=0,header=T)
#check plot iso
obj<-plotIso(obj,type="t5")

#check diff exp: this will become a DESeq2 obj,
#so any function can be apply to this
dds<-deIso(obj,formula=~condition,ref=T,iso5=T)
#plotMA
library(DESeq2)
plotMA(dds)
head(counts(dds))


#get counts
counts<-makeCounts(obj)
pls.obj<-isoPLSDA(obj)
isoPLSDAplot(pls.obj$components,obj@design[,1])
