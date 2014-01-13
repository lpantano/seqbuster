d<-getwd()
name<-paste(sep="",d,"/","temp/usRNA")
d<-read.table(name)

png("temp/ntcomposition.png",width=100,height=100,units="mm",res=150)

par(mfrow=c(2,2),cex=.5)

maxl<-max(d$V4)
minl<-min(d$V4)

l<-d

ma<-matrix(ncol=1,nrow=maxl-minl+1)
for (i in minl:maxl){
  #ma[i-17,1]<-length(l$V4[l$V2==i])
  ma[i-minl-1,1]<-sum(l$V3[l$V4==i])
}	
#print(ma)
barplot(t(ma[,1]),beside=T,col="blue4",xlab="size",ylab="expression",names.arg=minl:maxl,main=paste("Length distribution"))	 

l<-d[d$V1==1 ,]

ma<-matrix(ncol=4,nrow=maxl-minl+1)
for (i in minl:maxl){
  ma[i-minl-1,1]<-sum(l$V3[l$V4==i & l$V2=="A"])
  ma[i-minl-1,2]<-sum(l$V3[l$V4==i & l$V2=="T"])
  ma[i-minl-1,3]<-sum(l$V3[l$V4==i & l$V2=="G"])
  ma[i-minl-1,4]<-sum(l$V3[l$V4==i & l$V2=="C"])
  total<-sum(ma[i-minl-1,])
  ma[i-minl-1,]<-ma[i-minl-1,]/total

}	
#print(ma)

y=max(ma,na.rm=TRUE)
plot(minl:maxl,ma[,1],col="red",type="l",ylim=c(0,y+0.1*y),ylab="score",xlab="size",main=paste("First nt"))
points(minl:maxl,ma[,2],col="blue",type="l") 
points(minl:maxl,ma[,3],col="green",type="l") 
points(minl:maxl,ma[,4],col="yellow",type="l") 
legend("topleft",fill=c("red","blue","green","yellow"),legend=c("A","U","G","C"),border="white",horiz=TRUE,bty="n",cex=.7,title="size")

#midle
l<-d

ma<-matrix(ncol=4,nrow=maxl)
for (i in 1:maxl){
  ma[i,1]<-sum(l$V3[l$V1==i & l$V2=="A"])
  ma[i,2]<-sum(l$V3[l$V1==i & l$V2=="T"])
  ma[i,3]<-sum(l$V3[l$V1==i & l$V2=="G"])
  ma[i,4]<-sum(l$V3[l$V1==i & l$V2=="C"])
  total<-sum(ma[i,])
  ma[i,]<-ma[i,]/total
}	
#print(ma)

y=max(ma,na.rm=TRUE)
plot(1:maxl,ma[,1],col="red",type="l",ylim=c(0,y+0.1*y),ylab="score",xlab="position",main=paste("Along sequence"))
points(1:maxl,ma[,2],col="blue",type="l") 
points(1:maxl,ma[,3],col="green",type="l") 
points(1:maxl,ma[,4],col="yellow",type="l") 



dev.off()