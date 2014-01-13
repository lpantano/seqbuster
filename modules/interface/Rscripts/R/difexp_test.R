
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
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	#print(s)
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `id`,`seq`,`chr`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`DB`,`freq`,`amb` from `",s,"` ",options,"   `DB` like 'miRNA' AND `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	print(paste("miR-101 before norm",sum(temp$freq[temp$chr=="hsa-miR-101"])))
	print(temp[temp$chr=="hsa-miR-101" & temp$trimmed5==0 & temp$trimmed3==0 & temp$mut==0 & temp$addition3==0,])
	print(temp[1:2,])
	tempn<-applynorm(temp[,c(1,9)],list)
#print(temp[1:10,])
	temp<-temp[temp$amb<=3,1:9]
	print(paste("miR-101 after amb",sum(temp$freq[temp$chr=="hsa-miR-101"])))
    print(temp[temp$chr=="hsa-miR-101" & temp$trimmed5==0 & temp$trimmed3==0 & temp$mut==0 & temp$addition3==0,])
# 	print(tempn[1:10,])
	tempn<-merge(temp,tempn,by=1,all=FALSE)
		
	tempn<-tempn[tempn$DB==list$DB,]
#  	print(tempn[1:10,])
	temp<-tempn[,c(1,2,3,4,5,6,7,8,10)]
	names(temp)[9]<-"freq"
	print(paste("miR-101 after norm",sum(temp$freq[temp$chr=="hsa-miR-101"])))
  	print(temp[temp$chr=="hsa-miR-101" & temp$trimmed5==0 & temp$trimmed3==0 & temp$mut==0 & temp$addition3==0,])

	temp$freq<-round(temp$freq)
#  	print(temp[1:10,])
	query<-paste(sep="","select SUM(freq) AS freq from `",s,"` where `DB` like 'miRNA';")
	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
    print(max1[ns])
	table<-data.frame(id=0,freq=0)
	listmitemp<-unlist(unique(temp$chr))
	#listmitemp<-"hsa-miR-101"
	for (mi in listmitemp){
	
		pertemp<-temp[temp$chr==mi & temp$trimmed5==0 & temp$trimmed3==0 & temp$mut==0 & temp$addition3==0,]
		freqperfect<-0
		#print(pertemp)
		if (length(pertemp$chr)>0){freqperfect<-pertemp$freq}
		vartemp<-temp[temp$chr==mi & (temp$trimmed5!=0 | temp$trimmed3!=0 | temp$mut!=0 | temp$addition3!=0),]
		tempvalue<-sum(vartemp$freq)
		
		minvalue<-min(freqperfect,tempvalue,na.rm=T)
		if (mi=="hsa-miR-101"){
		  print(nrow(vartemp))
		  print (paste("minvalue: ",minvalue))
		}
		vartemp<-filter(vartemp,minvalue*cof,7,type="filter",freq=freqperfect,error=list$error)
		if (mi=="hsa-miR-101"){
		  print(nrow(vartemp))
		}
		
		if (length(vartemp$id)>0){

				parsetemp<-mapply(getid,vartemp$chr,vartemp$trimmed5,vartemp$trimmed3,vartemp$addition3,vartemp$mut)
				#devolver chrmutt5t3a3  and freq en una tabla add new values cbin()
# 				print (cutofftemp)
				freqtemp<-vartemp$freq
				if (length(pertemp$chr)>0)
				{
	# 			perfreqtemp<-0
					if (is.null(list$ref)==FALSE){
						labelper<-paste(sep="",pertemp$chr,",Ref")
					}else{
						labelper<-paste(sep="",pertemp$chr,",")
					}
					parsetemp<-append(as.character(parsetemp),as.character(paste(sep="",labelper)))
					freqtemp<-c(vartemp$freq,pertemp$freq)
# 					print (pertemp)
				}
# 				print (c(parsetemp,freqtemp))
				tabletemp<-data.frame(id=parsetemp,freq=freqtemp)
				
				listidtemp<-unlist(unique(parsetemp))
				freqtemp<-as.numeric(mapply(getfreq,listidtemp))
				table<-rbind(table,data.frame(id=listidtemp,freq=freqtemp))
# 				print (listidtemp)
		}else if(freqperfect>0 & is.null(list$ref)==F){
			labelper<-paste(sep="",mi,",Ref")
			table<-rbind(table,data.frame(id=labelper,freq=freqperfect))
		}
		
	}
# 	print (nrow(table))
#	listsamples1[[ns]]<-temp
# 	table$n<-round(table[,2])
	table[is.na(table)]<-1
	if (ns==1){
		all<-table[table$id!=0,]
		names(all)<-c('id',list$group1[ns])
		#tableall<-temp
		
	}else{
		
		names(table)<-c('id',list$group1[ns])
		all<-merge(all,table[table$id!=0,],by="id",all=TRUE)
		#tableall<-merge(tableall,temp,all=TRUE,by="id")
	}
}

#calculate pvalue intra groups
all[is.na(all)]<-1
table1<-all
table1all<-table1
names(table1all)<-c('id','t1')

#print(table1[grep("hsa-miR-101",table1$id) ,])

#print(paste("miR-101 after filtering",sum(table1[grep("hsa-miR-101",table1$id),2])))
#print (table1[1:10,])
#q()
cutoffsamples1<-(length(samplesnames1)-1)*outl;

tablepvalue<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
tableratio<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
if (length(list$group1)>1){
	colpvalue<-1
# 	listloci<-vector("list",length=nrow(table))
	
	for (s1 in 1:(length(list$group1)-1)){
		for (s2 in (s1+1):length(list$group1)){
			ratio<-max1[s1]/max1[s2]
			tablepvalue[,colpvalue]<-mapply(normtest,((all[,s1+1])),((all[,s2+1])),ratio)
			tableratio[,colpvalue]<-(all[,s1+1])/(all[,s2+1])
			colpvalue<-colpvalue+1
		}
	}
	#do coherent conditions
	if (colpvalue>1){
		print(outl)
		print(cutoffsamples1)
		tempratio<-apply(tableratio,1,median)
		temppvalue<-apply(tablepvalue,1,median)
		numpval<-apply(tablepvalue,1,numpvalfunc)
		numexp<-apply(table1[2:(length(samplesnames1)+1)],1,numexpfunc)
		
		table1$tempr<-tempratio
		table1$tempp<-temppvalue
		table1$tempnp<-numpval
		table1$tempnexp<-numexp
# 		print(table1[1:10,])
# 		table1<-all[(table1[,2]>=cs & table1[,2]<=cb) | table1[,3]>=pval,]
		table1all<-table1[,c(1,ncol(table1))]
		names(table1all)<-c('id','t1')
		write.table(table1,"temp/group1",sep="\t",quote=F,row.names=F)
 		table1<-table1[table1$tempnp<=cutoffsamples1 & table1$tempnexp<=outl,]
		#print(table1[1:10,])
				table1[,2]<-round(apply(table1[,2:(length(samplesnames1)+1)],1,median))
# 		print(table1[1:10,])
		table1<-table1[,1:2]
		
	}

}else{

	#temp<-apply(table[,2:(length(list$group1)+1)],1,mean)
	table1<-all
}
names(table1)<-c("id","freq1")


# table1[1:10,]

ns<-0
samplesnames2<-(list$group2)
listsamples2<-vector("list",length=length(list$group2))
max2<-1:length(list$group2)
#Load samples inside groups
for (s in list$group2){
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	#print(s)
	infogroup<-append(infogroup,paste(sep="","G2:",s))

	query<-paste(sep="","select `id`,`seq`,`chr`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`DB`,`freq`,`amb`  from `",s,"` ",options," `DB` like 'miRNA' AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	
	temp <- as.data.frame(fetch(rs))
	#print (temp[1:2,])
	tempn<-applynorm(temp[,c(1,9)],list)
	temp<-temp[temp$amb==1,1:9]
	tempn<-merge(temp,tempn,by=1,all=FALSE)
	tempn<-tempn[tempn$DB==list$DB,]
# 	print(tempn)
	
	temp<-tempn[,c(1,2,3,4,5,6,7,8,10)]
	names(tempn)[9]<-"freq"
	temp$freq<-round(temp$freq)
 	
	query<-paste(sep="","select SUM(freq) AS freq from `",s,"` where `DB` like 'miRNA';")
	print(query)
	rs <- dbSendQuery(con,query) 
	max2[ns]<- unlist(as.vector(fetch(rs)))
	listmitemp<-unlist(unique(temp$chr))
	table<-data.frame(id=0,freq=0)

	for (mi in listmitemp){
		pertemp<-temp[temp$chr==mi & temp$trimmed5==0 & temp$trimmed3==0 & temp$mut==0 & temp$addition3==0,]
		freqperfect<-0
		if (length(pertemp$chr)>0){freqperfect<-pertemp$freq}
		vartemp<-temp[temp$chr==mi & (temp$trimmed5!=0 | temp$trimmed3!=0 | temp$mut!=0 | temp$addition3!=0),]
		tempvalue<-sum(vartemp$freq)
		
		minvalue<-min(freqperfect,tempvalue,na.rm=T)
# 		print(vartemp)
		vartemp<-filter(vartemp,minvalue*cof,7,type="filter",freq=freqperfect,error=list$error)
# 		print(vartemp)
		
		if (length(vartemp$id)>0){

				parsetemp<-mapply(getid,vartemp$chr,vartemp$trimmed5,vartemp$trimmed3,vartemp$addition3,vartemp$mut)
				#devolver chrmutt5t3a3  and freq en una tabla add new values cbin()
				#print (parsetemp)
				freqtemp<-vartemp$freq
				if (length(pertemp$chr)>0)
				{
	# 			perfreqtemp<-0
					if (is.null(list$ref)==FALSE){
						labelper<-paste(sep="",pertemp$chr,",Ref")
					}else{
						labelper<-paste(sep="",pertemp$chr,",")
					}
					parsetemp<-append(as.character(parsetemp),as.character(paste(sep="",labelper)))
					freqtemp<-c(vartemp$freq,pertemp$freq)
				}
# 				print (parsetemp)
# 				print(freqtemp)
				tabletemp<-data.frame(id=parsetemp,freq=freqtemp)

				listidtemp<-unlist(unique(parsetemp))
				freqtemp<-as.numeric(mapply(getfreq,listidtemp))

				table<-rbind(table,data.frame(id=listidtemp,freq=freqtemp))
		}else if(freqperfect>0 & is.null(list$ref)==F){
			labelper<-paste(sep="",mi,",Ref")
			table<-rbind(table,data.frame(id=labelper,freq=freqperfect))
		}
		
		
		
	}

	table[is.na(table)]<-1
	if (ns==1){
		all<-table[table$id!=0,]
		names(all)<-c('id',list$group2[ns])
		#tableall<-temp
		
	}else{
		
		names(table)<-c('id',list$group2[ns])
		all<-merge(all,table[table$id!=0,],by="id",all=TRUE)
		#tableall<-merge(tableall,temp,all=TRUE,by="id")
	}
}

#calculate pvalue intra groups
all[is.na(all)]<-1
table2<-all
table2all<-table2
names(table2all)<-c('id','t2')

cutoffsamples2<-(length(samplesnames2)-1)*outl;


tablepvalue<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
tableratio<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
if (length(list$group2)>1){
	colpvalue<-1
# 	listloci<-vector("list",length=nrow(table))
	
	for (s1 in 1:(length(list$group2)-1)){
		for (s2 in (s1+1):length(list$group2)){
			ratio<-max2[s1]/max2[s2]
			tablepvalue[,colpvalue]<-mapply(normtest,((all[,s1+1])),((all[,s2+1])),ratio)
			tableratio[,colpvalue]<-(all[,s1+1])/(all[,s2+1])
			colpvalue<-colpvalue+1
		}
	}
	
	#do coherent conditions
	if (colpvalue>1){
		
		tempratio<-apply(tableratio,1,median)
		temppvalue<-apply(tablepvalue,1,median)
		numpval<-apply(tablepvalue,1,numpvalfunc)
		numpexp<-apply(table2[2:(length(samplesnames2)+1)],1,numexpfunc)
		table2$tempr<-tempratio
		table2$tempp<-temppvalue
		table2$tempnp<-numpval
		table2$tempnexp<-numpexp
# 		print(table1[1:10,])
# 		table2<-all[(table2[,2]>=cs & table2[,2]<=cb) | table2[,3]>=pval,]
		table2all<-table2[,c(1,ncol(table2))]
		names(table2all)<-c('id','t2')
		write.table(table2,"temp/group2",sep="\t",quote=F,row.names=F)
#		print(table2[1:10,])
		#table2<-table2[(table2$tempr>=cs & table2$tempr<=cb) & table2$tempnp>=cutoffsamples2 ,]
		table2<-table2[table2$tempnp<=cutoffsamples2 & table2$tempnexp<=outl,]
		
		table2[,2]<-round(apply(table2[,2:(ncol(table2)-2)],1,median))
# 		print(table1[1:10,])
		table2<-table2[,1:2]
		
	}
}else{

	#temp<-apply(table[,2:(length(list$group1)+1)],1,mean)
	table2<-all
}
names(table2)<-c("id","freq2")
###########remove not in sref
print ("after joing samples in group2")

#pvalue of the two groups
table2[,3]<-table2[,2]
table1[,3]<-table1[,2]
names(table2)<-c("id","freq2","freq22")
names(table1)<-c("id","freq1","freq11")

# print (table1[1:10,])
# print (table2[1:10,])
table<-merge(table1,table2,by="id",all=TRUE)

table[is.na(table)]<-0
checktable<-merge(table,table1all,by='id',all=TRUE)
checktable<-merge(checktable,table2all,by='id',all=TRUE)
print (checktable[1:10,])
checktable[is.na(checktable)]<-10-11	
checktable2<-checktable[(checktable[,2]==0 & checktable$t1<0 & checktable[,4]>0) | (checktable[,4]==0 & checktable$t2<0 & checktable[,2]>0) | (checktable[,2]>0 & checktable[,4]>0),]
#print (checktable2[1:10,])
table<-checktable2[,1:5]

ratio<-mean(max1)/mean(max2)
print(table[1:10,])

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
print ("after ztest")

if (list$pcorrect!="na"){
	sort<-sort(table$p,index.return=T)
	ind<-1:nrow(table)
	order<-unlist(lapply(ind,function (x) ind[sort$ix==x]))
	table$q<-mapply(BHcorrection,table$p,order,nrow(table))
}else{
	table$q<-1
}

# table[1:10,]
list$group2
list$group1
table$ratio<-table[,3]/table[,5]
# table[is.infinite(table)==T]
table$ratio[is.infinite(table$ratio)==T]<-9999
table$ratio[is.na(table$ratio)]<-1


table<-table[table[,3]!=0 | table[,5]!=0,]
print (table[1:10,])
tablename<-"difexp"
query<-paste(sep="","DROP TABLE IF EXISTS `",list$project,"`.`",tablename,"`;")
rs <- dbSendQuery(con,query) 

createtable<-paste(sep="","CREATE TABLE `",list$project,"`.`",tablename,"` (`chr` VARCHAR(36),`type` VARCHAR(36),`normG1` DOUBLE(8,2),`normG2` DOUBLE(8,2),`ratio` DOUBLE(8,2),pvalue DOUBLE(8,2),qvalue DOUBLE(8,2),INDEX (chr));")
rs <- dbSendQuery(con,createtable) 


for (r in 1:nrow(table)){
 	
	locus<-unlist(strsplit(table[r,1],","))[1]
	type<-unlist(strsplit(table[r,1],","))[2]
	loadtable<-paste(sep="","INSERT INTO `",list$project,"`.`",tablename,"` VALUES ('",locus,"','",type,"',",table[r,3],",",table[r,5],",",table$ratio[r],",",table$p[r],",",table$q[r],");")
#  	print (loadtable)
	rs <- dbSendQuery(con,loadtable) 
	
}

dbDisconnect(con) 

