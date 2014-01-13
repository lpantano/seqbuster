library("XML")
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")
list<-parse("temp/script_param")


MySQL(max.con = 1, fetch.default.rec = 10000000, force.reload = FALSE)
m <- dbDriver("MySQL")
con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw)

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
		id<-paste(sep="",id,"T5:",t5)
	}
	if (is.null(list$trimmed3)==FALSE & t3!=0  ){
		id<-paste(sep="",id,"T3:",t3)
	}
	if (is.null(list$addition3)==FALSE & a3!=0){
		id<-paste(sep="",id,"A3:",a3)
	}
	if (is.null(list$mut)==FALSE & m!=0){
	if (m!=0){
		
# 		mutation<-unlist(strsplit(gsub("[0-9]+","",m,""),""))
# 		nt1<-mutation[1]
# 		nt2<-mutation[2]
		pos<-as.numeric(gsub("[ATGC]+","",m))
		if (pos>=as.numeric(list$start) & pos<=as.numeric(list$end)){
# 			print (m)
			id<-paste(sep="",id,"M:",m)
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
# if (is.null(list$DB)==F){
# 	
# 	options<-paste(sep="",options," `DB` LIKE '",list$DB,"' AND")
# }
# if (is.null(list$locus)==F){
# 	
# 	options<-paste(sep="",options," `name` LIKE '%",list$locus,"%' AND")
# }
# if (is.null(list$ref)==F){
# 	
# 	options<-paste(sep="",options," `id` IN (select `id` from `",list$ref,"`) AND")
# }
cs<-0.5
cs<-as.numeric(list$cs)
cb<-1.5
cb<-as.numeric(list$cb)
pval<-0.05
pval<-as.numeric(list$pval)

temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])

cof<-as.numeric(list$cof)/100
# cof*2

#all chr of one sample
# list$reaf<-"CCF1"
# query<-paste(sep="","select `chr` AS freq from `",list$ref,"` ",options," group by `chr`;")
# print(query)
# rs <- dbSendQuery(con,query) 
# loci<- as.vector(fetch(rs))
# list$loci<-loci
infogroup<-vector()
scale<-as.numeric(list$norm)
ns<-0
listsamples1<-vector("list",length=length(list$group1))
max1<-1:length(list$group1)
table<-data.frame(id=0,freq=0)
listseqnames<-data.frame(seq='0',name='0')
#Load samples inside groups
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	#print(s)
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `seq`,`freq`,`name`,`DB` from `",s,"` ",options," ;")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",s,"` ",options," ;")
	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
# 	if (is.null(list$ref)==F){
# 		temp<-getcommon(temp)
# 	}
# 		print(temp[1:10,])
		table<-applynorm(temp[,c(1,2)],list)
# 		
		table<-merge(temp,table,by=1,all=FALSE)
		table<-table[table$DB==list$DB,1:3]
if (ns==1){
		all<-table[table$seq!=0,1:2]
		names(all)<-c('id',list$group1[ns])
		listseqnames<-temp[,c(1,3)]
		#tableall<-temp
		
	}else{
		
		names(table)<-c('id',list$group1[ns])
		all<-merge(all,table[table$seq!=0,1:2],by="id",all=TRUE)
		listseqnames<-merge(listseqnames,temp[,c(1,3)],by=1,all=TRUE)
# 		templistnames<-
		#tableall<-merge(tableall,temp,all=TRUE,by="id")
	}
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1
table1<-all

# print ((table1[1:10,]))
tablepvalue<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
tableratio<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
if (length(list$group1)>1){
	colpvalue<-1
# 	listloci<-vector("list",length=nrow(table))
	
	for (s1 in 1:(length(list$group1)-1)){
		for (s2 in (s1+1):length(list$group1)){
			ratio<-max1[s1]/max1[s2]
			tablepvalue[,colpvalue]<-mapply(statistic,((all[,s1+1])),((all[,s2+1])),ratio)
			tableratio[,colpvalue]<-(all[,s1+1])/(all[,s2+1])
			colpvalue<-colpvalue+1
		}
	}
	#do coherent conditions
	if (colpvalue>1){
		
		tempratio<-apply(tableratio,1,median)
		temppvalue<-apply(tablepvalue,1,median)
		table1[,2]<-tempratio
		table1[,3]<-temppvalue
# 		print(table1[1:10,])
# 		table1<-all[(table1[,2]>=cs & table1[,2]<=cb) | table1[,3]>=pval,]
		
		table1[,2]<-round(apply(all[,2:ncol(table1)],1,sum))
# 		print(table1[1:10,])
		table1<-table1[,1:2]
		
	}
# 	else{
# 		tempma<-tableratio
# 		tempma[,1]<-apply(tableratio,1,func)
# 		table1<-table[tempma[,1]==1,1:2]
# 		table1[,1]<-table[tempma[,1]==1,1]
# 		table1[,2]<-round(apply(table[tempma[,1]==1,c(2,3)],1,mean))
# 	}
}else{

	#temp<-apply(table[,2:(length(list$group1)+1)],1,mean)
	table1<-all
}
names(table1)<-c("id","freq1")
# table1[1:10,]
# print (nrow(table1))
# scale
#Load samples inside groups
# q()
ns<-0
listsamples2<-vector("list",length=length(list$group2))
max2<-1:length(list$group2)
#Load samples inside groups
for (s in list$group2){
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	#print(s)
	infogroup<-append(infogroup,paste(sep="","G2:",s))

	query<-paste(sep="","select `seq`,`freq`,`name`,`DB` from `",s,"` ",options," ;")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",s,"` ",options," ;")
	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
# 	if (is.null(list$ref)==F){
# 		temp<-getcommon(temp)
# 	}
# 		print(temp[1:10,])
		table<-applynorm(temp[,c(1,2)],list)
# 		
		table<-merge(temp,table,by=1,all=FALSE)
		table<-table[table$DB==list$DB,1:3]
# 	print(table[1:10,])
# 	print(listseqnames[1:2,])
	listseqnames<-merge(listseqnames,temp[,c(1,3)],by=1,all=TRUE)
		
	if (ns==1){
		all<-table[table$seq!=0,1:2]
		names(all)<-c('id',list$group2[ns])
		#tableall<-temp
		
	}else{
		
		names(all)<-c('id',list$group2[ns])
		all<-merge(all,table[table$seq!=0,1:2],by="id",all=TRUE)
		#tableall<-merge(tableall,temp,all=TRUE,by="id")
	}
}

#calculate pvalue intra groups
all[is.na(all)]<-0
table2<-all

# print (nrow(table2))
tablepvalue<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
tableratio<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
if (length(list$group2)>1){
	colpvalue<-1
# 	listloci<-vector("list",length=nrow(table))
	
	for (s1 in 1:(length(list$group2)-1)){
		for (s2 in (s1+1):length(list$group2)){
			ratio<-max2[s1]/max2[s2]
			tablepvalue[,colpvalue]<-mapply(statistic,((all[,s1+1])),((all[,s2+1])),ratio)
			tableratio[,colpvalue]<-mapply(norm,((all[,s1+1])),((all[,s2+1])),max2[s1],max2[s2])
			colpvalue<-colpvalue+1
		}
	}
	
	#do coherent conditions
	if (colpvalue>1){
		
		tempratio<-apply(tableratio,1,median)
		temppvalue<-apply(tablepvalue,1,median)
		table2[,2]<-tempratio
		table2[,3]<-temppvalue
# 		print(table1[1:10,])
# 		table2<-all[(table2[,2]>=cs & table2[,2]<=cb) | table2[,3]>=pval,]
		
		table2[,2]<-round(apply(all[,2:ncol(table2)],1,sum))
# 		print(table1[1:10,])
		table2<-table2[,1:2]
		
	}
}else{

	#temp<-apply(table[,2:(length(list$group1)+1)],1,mean)
	table2<-all
}

names(table2)<-c("id","freq2")
###########remove not in sref
# print (table1[1:10,])
# print (table2[1:10,])
#pvalue of the two groups
table1[,3]<-table1[,2]
table2[,3]<-table2[,2]
table<-merge(table1,table2,by='id',all=TRUE)

table[is.na(table)]<-0
# print (nrow(table))
# if (is.null(list$sref)==F){
# 	table<-table[table[,lsref]>0,]
# 
# }
# table

ratio<-mean(max1)/mean(max2)
# print (table[1:10,])
if (list$st=="bayesian"){
	table$p<-mapply(statistic,table[,3],table[,5],ratio)
}else if (list$st=="binomial"){
	table$p<-mapply(binomial,table[,3],table[,5],mean(max1),mean(max2))
}else if (list$st=="ztest"){
	table$p<-mapply(ztest,table[,3],table[,5],sum(table[,3]),sum(table[,5]))
}else if (list$st=="fisher"){
	table$p<-mapply(fisher,table[,3],table[,5],sum(table[,3]),sum(table[,5]))
}else{
	table$p<-1
}
table$q<--1

if (list$pcorrect!="na"){
	sort<-sort(table$p,index.return=T)
	ind<-1:nrow(table)
	order<-unlist(lapply(ind,function (x) ind[sort$ix==x]))
	table$q<-mapply(BHcorrection,table$p,order,nrow(table))
}else{
	table$q<-1
}
# print (order)
# print (sort)
# print(table[1:10,])
list$group2
list$group1
table$ratio<-table[,3]/table[,5]
# table[is.infinite(table)==T]
table$ratio[is.infinite(table$ratio)==T]<-9999
table$ratio[is.na(table$ratio)]<-1
# q()
des<-nrow(table[table$q<0.01 & (table$ratio>=1.5 | table$ratio<=0.5),])
nrow(table)
# des
table<-table[table[,3]!=0 & table[,4]!=0,]
# max2
#  ip<-unlist(strsplit(args[1],"ip"))
tablename<-"result"

query<-paste(sep="","DROP TABLE IF EXISTS `",list$project,".",tablename,"`;")
rs <- dbSendQuery(con,query) 

createtable<-paste(sep="","CREATE TABLE `",list$project,".",tablename,"` (`chr` VARCHAR(36),`type` VARCHAR(36),`normG1` DOUBLE(8,2),`normG2` DOUBLE(8,2),`ratio` DOUBLE(8,2),pvalue DOUBLE(8,2),qvalue DOUBLE(8,2),INDEX (chr));")
rs <- dbSendQuery(con,createtable) 


for (r in 1:nrow(table)){

	
	nameseq<-listseqnames[listseqnames$seq==table[r,1],2:ncol(listseqnames)]
	nameseq<-nameseq[nameseq!="0"]
# 	print(nameseq[1])
	locus<-nameseq[1]
	type<-table[r,1]
	loadtable<-paste(sep="","INSERT INTO `",list$project,".",tablename,"` VALUES ('",locus,"','",type,"',",table[r,3],",",table[r,5],",",table$ratio[r],",",table$p[r],",",table$q[r],");")
#  	print (loadtable)
	rs <- dbSendQuery(con,loadtable) 
	
}


dbDisconnect(con) 
