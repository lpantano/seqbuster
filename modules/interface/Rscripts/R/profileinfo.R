d<-getwd()

g1<-read.table(paste(sep="",d,"/","temp/group1"),header=T)
g2<-read.table(paste(sep="",d,"/","temp/group2"),header=T)
name<-paste(sep="",d,"/","temp/profile.info")
d<-read.table(name)
id<-d[1,1]
print(id)

gs1<-g1[g1[,1]==id,]
gs2<-g2[g2[,1]==id,]

l1<-ncol(g1)-4
l2<-ncol(g2)-4

x1<-1:(l1-1)
x2<-1:(l2-1)
y1<-log(as.numeric(gs1[1,2:l1])+1)
y2<-log(as.numeric(gs2[1,2:l2])+1)

maxylim<-max(y1,y2)+1
minylim<-min(y1,y2)-1


print(gs1)
print(gs2)
print(x1)

write.table(gs1,"temp/group1table",row.names=F,quote=F)
write.table(gs2,"temp/group2table",row.names=F,quote=F)



jpeg("temp/profile.jpg",res=100)
plot(x1,y1,type='l',ylab="expression (log)",xlab="",col="green",xlim=c(0,max(l1,l2)),ylim=c(minylim,maxylim))
points(x2,y2,type='l',ylab="expression",xlab="",col="blue")
legend("bottomright",c(paste(names(g1)[2:l1],collapse=","),paste(names(g2)[2:l2],collapse=",")),fill=c("green","blue"),cex=0.7)
dev.off()