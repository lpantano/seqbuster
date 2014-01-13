

g2<-read.table("temp/random_g2")
g1<-read.table("temp/random_g1")


png("temp/random-srna.png")

p<-round(unlist(t.test(g2[,1],g1[,1])[3]),digits=3)

h<-hist(g2[,1],col=rgb(0,0,1,1/4),border=NA,freq=F,main="Secondary structure distances",xlab="distance score")
hist(g1[,1],add=T,col=rgb(0,1,0,1/4),border=NA,freq=F)
legend("topleft",fill=c("blue4","green4"),legend=c("random structures group","sRNA and random structure group"),border="white",bty=0,cex=.7)
my<-max(h$density)
text(10,my/2,labels=paste("p-val= ",p))
dev.off()