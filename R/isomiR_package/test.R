files<-c("/home/lpantano/crickhome/isomirs/mirbase20/y0d204.hsa.fa.ad.mirna.out.freq.par",
         "/home/lpantano/crickhome/isomirs/mirbase20/y0d2.hsa.fa.ad.mirna.out.freq.par",
         "/home/lpantano/crickhome/isomirs/mirbase20/y0d34.hsa.fa.ad.mirna.out.freq.par",
         "/home/lpantano/crickhome/isomirs/mirbase20/y0d4.hsa.fa.ad.mirna.out.freq.par")

d<-data.frame(condition=c(1,1,2,2))
obj<-loadIso(files,d)
obj<-plotIso(obj,type="t5")
dds<-deIso(obj,formula=~condition)
plotTop(dds)

f<-read.table(files[1],skip=1,header=F)