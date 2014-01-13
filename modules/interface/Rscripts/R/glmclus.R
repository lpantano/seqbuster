library("XML")
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")



normtest<-function(x1,x2,n1,n2){
	
	p<-pnorm(x1,mean=x2,sd=x2/3)
	if(p>0.9){
		p<-1-p
	}
	return(p)
}
cs<-0.5
cb<-1.5
pval<-0.05

numratiofunc<-function(v){
	
	n<-length(v[v<cb & v>=cs])
	return(n)
}


numpvalfunc<-function(v){
	
	n<-length(v[v<=pval])
	return(n)
}

numexpfunc<-function(v){
	
	n<-length(v[v<=1])
	return(n)
}


list<-parse("temp/script_param")
outl<-as.numeric(list$outl)
list$DB="miRNA"

MySQL(max.con = 1, fetch.default.rec = 10000000, force.reload = FALSE)
m <- dbDriver("MySQL")
if (port==0){
	con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw)
}else{
	con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw,port=port)
}


norm<-function (n1,n2,max1,max2){
	
	r<-(n1/max1*1000000)/(n2/max2*1000000)
	return (r)
}


func<-function(n){
	p<-0
	if (n<=list$ratiocutoff[1] & n>=list$ratiocutoff[2]){
		p<-1
	}
	return (p)
}

normalization<-function(n,max){
	
	n<-as.numeric(n)
	nn<-round(n/max*scale)
	if (nn==0){nn=1}
	return(nn)
}

getid<-function(chr,t5,t3,a3,m){
 	id<-paste(sep="",chr,",")
	
	if (is.null(list$trimmed5)==FALSE & t5!=0  ) {
		id<-paste(sep="",id," tr5:",t5)
	}
	if (is.null(list$trimmed3)==FALSE & t3!=0  ){
		id<-paste(sep="",id," tr3:",t3)
	}
	if (is.null(list$addition3)==FALSE & a3!=0){
		id<-paste(sep="",id," ad3:",a3)
	}
	if (is.null(list$mut)==FALSE & m!=0){
		if (m!=0){
			
			
			pos<-as.numeric(gsub("[ATGC]+","",m))
			if (pos>=as.numeric(list$start) & pos<=as.numeric(list$end)){
# 			print (m)
				id<-paste(sep="",id," s:",m)
			}
			
		}
	}
	
	return(id)
	
}

getfreq<-function(id){
	
	freq<-sum(tabletemp$freq[tabletemp$id==id])
# 	print (tabletemp$id==id)
	return (freq)
	
}
if (is.null(list$freq1)==F){
	temp<-unlist(strsplit(list$freq1," "))
	options<-paste(sep=" ","where `freq` <=",list$freq1," AND")
}else{
	options<-"where `freq` > 0 AND"
}
if (is.null(list$freq2)==F){
	temp<-unlist(strsplit(list$freq1," "))
	options<-paste(sep=" "," `freq` >=",list$freq2," AND")
}
if (is.null(list$len1)==F){
	temp<-unlist(strsplit(list$len1," "))
	options<-paste(sep=" ","  `len` <=",list$len1," AND")
}
if (is.null(list$len2)==F){
	temp<-unlist(strsplit(list$len2," "))
	options<-paste(sep=" "," `len` >=",list$len2," AND")
}

if (is.null(list$locus)==F){
	
	options<-paste(sep="",options," `chr` LIKE '%",list$locus,"%' AND")
}

cs<-0.5
#cs<-as.numeric(list$cs)
cb<-1.5
#cb<-as.numeric(list$cb)
pval<-0.05
#pval<-as.numeric(list$pval)

numpvalfunc<-function(v){
	
	n<-length(v[v<=pval])
	return(n)
}

numexpfunc<-function(v){
	
	n<-length(v[v<=1])
	return(n)
}


temp<-unlist(strsplit(options," "))
if (length(temp)>=2){
	
	options<-paste(collapse=" ",temp[1:(length(temp)-1)])
	options<-paste(options,"AND")
	
}
cof<-as.numeric(list$cof)/100

infogroup<-vector()
scale<-as.numeric(list$norm)
ns<-0
listsamples1<-vector("list",length=length(list$group1))
samplesnames1<-(list$group1)
max1<-1:length(list$group1)
table<-data.frame(id=0,freq=0)
#Load samples inside groups
for (s in list$group1){
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

	table<-applynorm(temp[,c(1,2)],scaleval,q1val,q2val,typetrval,typetdval)
	table<-merge(temp,table,by=1,all=FALSE)
	
	if (ns==1){
		all<-table[,c(1,3)]
		names(all)<-c('id',[ns])
# 		print(all[1:10,])
		#tableall<-temp
		
	}else{
		table<-table[,c(1,3)]
		names(table)<-c('id',samplesnames1[ns])
		all<-merge(all,table,by="id",all=TRUE)
		
	}
	
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1
table1<-all
write.table(table1,"temp/group1",sep="\t",quote=F,row.names=F)


ns<-0
samplesnames2<-(list$group2)
listsamples2<-vector("list",length=length(list$group2))
max2<-1:length(list$group2)
#Load samples inside groups
for (s in list$group2){
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

		table<-applynorm(temp[,c(1,2)],scaleval,q1val,q2val,typetrval,typetdval)
# 		
		table<-merge(temp,table,by=1,all=FALSE)
			
	if (ns==1){
		all<-table[,c(1,3)]
		names(all)<-c('id',samplesnames2[ns])

	}else{
		table<-table[,c(1,3)]
		names(table)<-c('id',samplesnames2[ns])
		all<-merge(all,table,by="id",all=TRUE)
		
	}
}

#calculate pvalue intra groups
all[is.na(all)]<-1
table2<-all

write.table(table2,"temp/group2",sep="\t",quote=F,row.names=F)

l1<-length(max1)+1
l2<-length(max2)+1
all<-merge(table1[,1:l1],table2[,1:l2],by=1,all=TRUE)
all[is.na(all)]<-0

n<-length(max1)

for (i in 1:nrow(all)){
	v<-all[i,]
#	print (v)
	gs1<-v[2:l1]
	gs2<-v[(l1+1):(l2+l1-1)]
	len1<-length(gs1)
	len2<-length(gs2)
	n<-len1+len2
	
	table<-matrix(ncol=3,nrow=n)
	values<-unlist(c(gs1,gs2))
#print(values)
#	print(c(rep(0,len1),rep(1,len2)))
#	print(rep(1,n))
	table<-data.frame(case=c(rep(0,len1),rep(1,len2)),count=values,norm=rep(1,n))
	table$norm<-unlist(c(max1,max2))
	table$f<-table$norm/min(table$norm)
	de<-glm(count~as.factor(case)+offset(log(f)),data=table,family="quasipoisson")
	su<-summary(de)
	all$pval[i]<-su$coefficients[2,4]
}

libsize<-unlist(c(max1,max2))
for (i in 1:(length(libsize))){
	all[,i+1]<-all[,i+1]/libsize[i]*min(libsize)
	all[,i+1]<-round(all[,i+1],digits=2)
}
write.table(all,paste(list$path,"table.txt",sep=""),sep="\t",quote=F,row.names=F)

doc = newXMLDoc()
top<-newXMLNode("html")
br<-newXMLNode("p","Differential expression with GLM",attrs=c(align="center"),parent=top)
a<-newXMLNode("a","see table",attrs=c(title="difexp",href="table.txt",parent=br))
saveXML(top,file=paste(sep="",list$path,"result.html"))

dbDisconnect(con) 



#y<-matrix(ncol=14,nrow=nrow(all))

#y[,1:14]<-as.matrix(all[,2:15]


#d <- DGEList(counts=y,group=x,lib.size=c(lib1,lib2),genes=all[,1])

#d <- d[rowSums(d$counts) >= 5, ]
#rownames(dP$counts)<-as.character(all[,1])
#d <- estimateCommonDisp(d)
#d <- estimateTagwiseDisp(d)

#de<-exactTest(d,common.disp=FALSE)
#top<-topTags(de,n=nrow(all))$table

#temp<-merge(all,top,by=1,all=FALSE)

