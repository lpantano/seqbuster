
library("XML")
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")
list<-parse("temp/script_param")


MySQL(max.con = 1, fetch.default.rec = 10000000, force.reload = FALSE)
m <- dbDriver("MySQL")
if (port==0){
	con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw)
}else{
	con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw,port=port)
}


direchetli<-function(n,m,k,a){

	const<-lgamma(n+a)/lgamma(a)^2
	p<-const*lgamma(m+a)
	return(p)
}


direchetls<-function(n,m,k,a){

	const<-lgamma(n+a+m)
	p<-const/lgamma(a)
	return(p)
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
gettype<-function(name){

	split<-unlist(strsplit(name,","))
	if (length(split)>1){
		t<-1
# 		print (split)
	}else{
		t<-0
	}
	if (list$div==0){
		t<-1
	}
	return(t)
}
getid<-function(chr,t5,t3,a3,m){
 	id<-paste(sep="",chr,",")
	ok<-0
	if (is.null(list$trimmed5)==FALSE & t5!=0  ) {
		id<-paste(sep="",id,"T5:",t5)
		ok<-1
	}
	if (is.null(list$trimmed3)==FALSE & t3!=0  ){
		id<-paste(sep="",id,"T3:",t3)
		ok<-1
	}
	if (is.null(list$addition3)==FALSE & a3!=0){
		id<-paste(sep="",id,"A3:",a3)
		ok<-1
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
			ok<-1
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

if (length(list$group1)<3){
	error<-"3 or more samples are needed to run the analysis...!"
	
}


if (is.null(list$linker)==F){
	
	options<-paste(sep="",options," `tag` LIKE '",list$linker,"' AND")
}
if (is.null(list$DB)==F){
	
	options<-paste(sep="",options," `DB` LIKE '",list$DB,"' AND")
}
if (is.null(list$locus)==F){
	
	options<-paste(sep="",options," `chr` LIKE '%",list$locus,"%' AND")
}
# if (is.null(list$sref)==F){
# 	
# 	options<-paste(sep="",options," `id` IN (select `id` from `",list$sref,"`) AND")
# }
list$div<-0
if (is.null(list$trimmed5)==FALSE ) {
		list$div<-1
	}
if (is.null(list$trimmed3)==FALSE ){
		list$div<-1
	}
if (is.null(list$addition3)==FALSE){
		list$div<-1
	}
if (is.null(list$mut)==FALSE){
		list$div<-1
}
if (is.null(list$ref)==FALSE ){
		list$div<-1
}
alpha<-as.numeric(list$alpha)

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

scale<-as.numeric(list$norm)
ns<-0
listsamples1<-vector("list",length=length(list$group1))
max1<-1:length(list$group1)
table<-data.frame(id=0,freq=0)
#Load samples inside groups
for (s in list$group1){
	sname<-s
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1

	#print(s)
	query<-paste(sep="","select `id`,`seq`,`chr`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`freq` AS freq from `",s,"` ",options," AND `amb`=1 AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+' AND `mut` NOT REGEXP '^[0-9]+[ACTG]+[0-9]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	
# 	if (is.null(list$sref)==F){
# 		temp<-getcommon(temp,nloci[1,1])
# 	}
	
	tempn<-applynorm(temp[,c(1,8)],list)
#  	print(temp[1:10,])
	tempn<-merge(temp,tempn,by=1,all=FALSE)
# 	print(nrow(temp))
	temp<-tempn[,c(1,2,3,4,5,6,7,9)]
	names(temp)[8]<-"freq"
	temp$freq<-round(temp$freq)
#  	print(temp[1:10,])
# 	print(sum(temp$freq))
	query<-paste(sep="","select SUM(freq) AS freq from `",s,"` ",options," AND `amb`=1 group by `DB`;")
# 	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
	table<-data.frame(id=0,freq=0)
	listmitemp<-unlist(unique(temp$chr))
# 	print (listmitemp)
	for (mi in listmitemp){
	
		pertemp<-temp[temp$chr==mi & temp$trimmed5==0 & temp$trimmed3==0 & temp$mut==0 & temp$addition3==0,]
		freqperfect<-0
		if (length(pertemp$chr)>0){freqperfect<-pertemp$freq}
		vartemp<-temp[temp$chr==mi & (temp$trimmed5!=0 | temp$trimmed3!=0 | temp$mut!=0 | temp$addition3!=0),]
		tempvalue<-sum(vartemp$freq)
		
		minvalue<-min(freqperfect,tempvalue,na.rm=T)
		
		vartemp<-filter(vartemp,minvalue*cof,7,type="filter",freq=freqperfect,error=list$error)
# 		print(listmitemp)
		
		if (length(vartemp$id)>0){
# 			for (nvar in 1:length(vartemp$id)){
# 				totalfreqid<-sum(vartemp$freq)
				#get all freq variation and cutoff 10%	
# 				cutofftemp<-cof*totalfreqid
# 				vartemp<-vartemp[vartemp$freq>=cutofftemp,]
				parsetemp<-mapply(getid,vartemp$chr,vartemp$trimmed5,vartemp$trimmed3,vartemp$addition3,vartemp$mut)
				
				#devolver chrmutt5t3a3  and freq en una tabla add new values cbin()
# 				print (cutofftemp)
# 				parsetemp<-parsetemp[parsetemp!=0]
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
# 				if (mi=="hsa-miR-411"){
# 					print (pertemp$freq)
# 					print (listidtemp)
# 					print (freqtemp)
# # 					q()
# 				}
				type<-mapply(gettype,listidtemp)
				temptable<-data.frame(id=listidtemp,type=type,freq=freqtemp)
# 				print(temptable)
				temptable<-temptable[temptable$type!=0,c(1,3)]
# 				print(temptable)
				table<-rbind(table,temptable)
				
# 				print (listidtemp)
		}
		
		
		
	}
# 	print(table[1:10,])
# 	table<-applynorm(table,list)
# 	print (table[1:10,])
	
#	listsamples1[[ns]]<-temp
# 	table$freq<-round(table[,2])*100
	if (ns==1){
		all<-table[table$id!=0,]
		names(all)<-c('id',sname)
		#tableall<-temp
		
	}else{
		
		names(table)<-c('id',sname)
		all<-merge(all,table[table$id!=0,],by="id",all=TRUE)
		#tableall<-merge(tableall,temp,all=TRUE,by="id")
	}
# 	print (all[1:10,])
}
all[is.na(all)]<-0
# c<-apply(all[,2:ncol(all)],1,min)
# all<-all[c>0,]
table<-all


if (list$clustme!="standar"){
###############################################
ma<-matrix(ncol=ncol(table)-1,nrow=ncol(table)-1)
ind1<-0
# alpha
# table
for (n1 in 2:(ncol(table))){
	nmi<-table[table[,n1]>0,1]
	n<-length(nmi)
	ind1<-ind1+1
	ind2<-0
	for (n2 in (2):ncol(table)){
		ind2<-ind2+1
		mmi<-table[table[,n2]>0,1]
		m<-length(mmi)
		k<-length(unlist(unique(c(mmi,nmi))))
		mu<-lgamma(k*alpha)^2
		de<-lgamma(n+k*alpha)*lgamma(m+k*alpha)
		const<-mu/de
		values<-mapply(direchetli,table[,n1],table[,n2],alpha,k)
		Li<-const*prod(values)
		
		mu<-lgamma(k*alpha)
		de<-lgamma(n+m+k*alpha)
		const<-mu/de
		values<-mapply(direchetls,table[,n1],table[,n2],alpha,k)
		Ls<-const*prod(values)
		
		d<-log(Li/(Ls+Li))
		ma[ind2,ind1]<-d
	}
}
# ma
row.names(ma)<-list$group1
names(ma)<-list$group1
ma<-scale(ma)
# ma
d<-as.dist(ma)
clust<-hclust(d,method="complete")
jpeg(paste(sep="",list$path,"img.jpg"),quality = 100)
par(mar=c(1,1,1,1))
plot(clust)
dev.off()
###############################################################
}else{
#####################################################
d<-as.dist(1-cor((table[,2:ncol(table)]),method="spearman"))
# d
clust<-hclust(d,method="complete")
# clust
jpeg(paste(sep="",list$path,"img.jpg"),quality = 100)
par(mar=c(1,1,1,1))
plclust(clust,labels=list$group1)
dev.off()
######################################################
}

doc = newXMLDoc()
top=newXMLNode("div",attrs=c(style="text-align: center;"))
newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = top)
save<-saveXML(top,file=paste(sep="",list$path,"result.html"))

dbDisconnect(con) 
