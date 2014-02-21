files<-c("/home/lpantano/crickhome/isomirs/mirbase20/y0d204.hsa.fa.ad.mirna.out.freq.par",
         "/home/lpantano/crickhome/isomirs/mirbase20/y0d2.hsa.fa.ad.mirna.out.freq.par",
         "/home/lpantano/crickhome/isomirs/mirbase20/y0d34.hsa.fa.ad.mirna.out.freq.par",
         "/home/lpantano/crickhome/isomirs/mirbase20/y0d4.hsa.fa.ad.mirna.out.freq.par")
library(isomiRs)
d<-data.frame(condition=c("p","p","c","c"))
row.names(d)<-paste(d[,1],1:2,sep="")

obj<-loadIso2(files,d)
#check plot iso
obj<-plotIso(obj,type="sub")

#check diff exp
dds<-deIso(obj,formula=~condition,ref=T)
plotTop(dds)

f<-read.table(files[1],skip=1,header=F)

head(obj@expList[["p1"]])