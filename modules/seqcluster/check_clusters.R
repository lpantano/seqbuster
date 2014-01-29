#this script read files created by seqcluster
#setwd(output seqcluter)
#and run the script. 
#You will get how complex are your clusters
#and how reliable are your results

library(igraph)
tls<-read.table("nodes.txt",fill=TRUE,header=F)
c<-read.table("nodesInfo.sort")

g <- graph.data.frame(tls, directed = F, vertices = c)
plot(g, layout = layout.fruchterman.reingold, 
     vertex.color = V(g)$V2, vertex.size = 2, 
     vertex.label = NA, 
     edge.arrow.size = 0.3, vertex.frame.color = V(g)$V2)


tls.p<-tls[!duplicated(tls[,1:2]),]
counts.o<-as.data.frame(table(tls.p[,1]))
tls.merge<-merge(tls.p,counts.o,by=1)
tls.o<-tls.p[order(tls.p$V2,tls.p$V4,decreasing=T),]
tls.o<-tls.o[!duplicated(tls.o$V2),]
tls.merge<-merge(tls.merge,tls.o[,c(2,4)])
tls.merge$ratio<-tls.merge$V3/tls.merge$V4

tls.merge.o<-tls.merge[order(tls.merge$V1,tls.merge$ratio,decreasing=T),]

d<-density((tls.merge[tls.merge$Freq>1 & tls.merge$V3!=-1,"ratio"]))
plot(d, main="Importance of original clusters to new clusters",
     xlab="% of sequences given to new cluster",
     ylab="# of original clusters")
polygon(d, col="steelblue", border="blue") 


