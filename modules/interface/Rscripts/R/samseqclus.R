
library("XML")
library("RMySQL")
library("samr")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")

d<-getwd()
name<-paste(sep="",d,"/","temp/args.txt")
d<-read.table(name)
opt<-as.vector(d[,1])
samplesnames1<-unlist(strsplit(opt[2],":"))
samplesnames2<-unlist(strsplit(opt[3],":"))
projectval<-opt[4]
typetrval<-"log2"
typetval<-"scale"
scaleval<-1000000
path<-opt[1]

MySQL(max.con = 16, fetch.default.rec = 5000000, force.reload = TRUE)

m <- dbDriver("MySQL")

if (port==0){
	con <- dbConnect(m,host=hostname,user=username,db=projectval,password=pssw)
}else{
	con <- dbConnect(m,host=hostname,user=username,db=projectval,password=pssw,port=port)
}

infogroup<-vector()
ns<-0
listsamples1<-vector("list",length=length(samplesnames1))
max1<-1:length(samplesnames1)
table<-data.frame(id=0,freq=0)
#Load samples inside groups
for (s in samplesnames1){
	#s<-paste(sep="",s,projectval)
	ns<-ns+1
	print(s)
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `idu`,`freq` from `",projectval,"`.`",s,"clusmap` where`idu`>0 group by `idu`;")
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",projectval,"`.`",s,"`  ;")
# 	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))

	temp[,2]<-temp[,2]/max1[ns]*scaleval
	table<-temp
	
	if (ns==1){
		all<-table
		names(all)<-c('id',samplesnames1[ns])

		
	}else{
		table<-table
		names(table)<-c('id',samplesnames1[ns])
		all<-merge(all,table,by="id",all=TRUE)
		
	}
	
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1
table1<-all

ns<-0

listsamples2<-vector("list",length=length(samplesnames2))
max2<-1:length(samplesnames2)
#Load samples inside groups
for (s in samplesnames2){
	#s<-paste(sep="",s,projectval)
	ns<-ns+1
	print(s)
	infogroup<-append(infogroup,paste(sep="","G2:",s))

	query<-paste(sep="","select `idu`,`freq` from `",projectval,"`.`",s,"clusmap` where `idu`>0 group by `idu`;")
#  	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",projectval,"`.`",s,"` ;")
#  	print(query)
	rs <- dbSendQuery(con,query) 
	max2[ns]<- unlist(as.vector(fetch(rs)))
	temp[,2]<-temp[,2]/max1[ns]*scaleval
	table<-temp			
	if (ns==1){
		all<-table
		names(all)<-c('id',samplesnames2[ns])

	}else{
		table<-table
		names(table)<-c('id',samplesnames2[ns])
		all<-merge(all,table,by="id",all=TRUE)
		
	}
}

#calculate pvalue intra groups
all[is.na(all)]<-1
table2<-all

all<-merge(table1,table2,by=1,all=TRUE)
all[is.na(all)]<-0

l1<-length(max1)
l2<-length(max2)

x<-as.matrix(all[,2:ncol(all)])
y<-c(rep(1,l1),rep(2,l2))

data=list(x=x,y=y, geneid=as.character(all[,1]),genenames=paste("cluster-",as.character(all[,1]),sep=""), logged2=TRUE)
samr.obj<-samr(data, resp.type="Two class unpaired", nperms=100)

delta.table <- samr.compute.delta.table(samr.obj)
deltaminus<- delta.table[delta.table[,5]<=0.05,1]

#print(delta.table[1:10,])

deltaminus<-sort(deltaminus)
delta<-deltaminus[1]
print(delta)
samr.plot(samr.obj,delta)

siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, data, delta.table)

samr.assess.samplesize.obj<- samr.assess.samplesize(samr.obj, data, log2(1.5))

png(paste(path,"/",projectval,"/samplesize.png",sep=""))
samr.assess.samplesize.plot(samr.assess.samplesize.obj)
dev.off()

write.table(siggenes.table$genes.up,paste(path,"/",projectval,"/tableup.txt",sep=""),sep="\t",quote=F,row.names=F)

write.table(siggenes.table$genes.lo,paste(path,"/",projectval,"/tabledown.txt",sep=""),sep="\t",quote=F,row.names=F)

d<-getwd()
write.table(table,paste(sep="",d,"/","temp/df.done"))
dbDisconnect(con) 


