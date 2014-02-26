
library(isomiRs)


#setwd("path2seqbuster/R/isomiR_package/test")

#seq path to source code
files<-c("y0d2.hsa.fa.ad.new.mirna",
          "y0d34.hsa.fa.ad.new.mirna",
         "y66d0.hsa.fa.ad.new.mirna",
         "y80d0.hsa.fa.ad.new.mirna"
         )


d<-data.frame(condition=c("p","p","c","c"))
row.names(d)<-paste(d[,1],1:2,sep="")

obj<-loadIso2(files,d,skip=0,header=T)
#check plot iso
obj<-plotIso(obj,type="t5")

#check diff exp: this will become a DESeq2 obj,
#so any function can be apply to this
dds<-deIso(obj,formula=~condition,ref=T,iso5=T)
#plotMA
library(DESeq2)
plotMA(dds)
head(counts(dds))
