d<-getwd()
name<-paste(sep="",d,"/","temp/cluster.info")
d<-read.table(name)

if (nrow(d)>1){
	

fin<-1
temppos<-d$V2[1]
for (i in 2:nrow(d)){
	if (temppos < d$V1[i]){
		fin<-i-1
	}
	
	temppos<-d$V2[i]
}
	d<-d[1:fin,]
}

print(d)

minp<-min(d$V1)
maxp<-max(d$V2)
p<-seq(minp,maxp,1)
table<-data.frame(pos=p,exp=rep(0,length(p)))
for (i in p){
	exp<-sum(d$V3[d$V1<=i & d$V2>=i])
	table$exp[table$pos==i]<-exp
}
label<-1:length(p)
table$p<-label
jpeg("temp/cluster.jpg")
plot(table$p,table$exp,type='l',ylab="expression",xlab="position")
dev.off()