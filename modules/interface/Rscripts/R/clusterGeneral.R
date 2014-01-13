
library("XML")
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")
list<-parse("temp/script_param")


MySQL(max.con = 1, fetch.default.rec = 10000000, force.reload = FALSE)
m <- dbDriver("MySQL")
con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw)




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
}if (is.null(list$linker)==F){
	
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

#all chr of one sample
# list$reaf<-"CCF1"
# query<-paste(sep="","select min(freq) as freq AS freq from `",list$ref,"` ",options," group by `DB`;")
# print(query)
# rs <- dbSendQuery(con,query) 
# minfreq<- as.vector(fetch(rs))

scale<-as.numeric(list$norm)
ns<-0
listsamples1<-vector("list",length=length(list$group1))
max1<-1:length(list$group1)

#Load samples inside groups

table<-data.frame(id=0,freq=0)
#Load samples inside groups
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	#print(s)

	query<-paste(sep="","select `seq`,`freq` from `",s,"` ",options," AND `amb`=1 ;")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",s,"` ",options," AND `amb`=1 group by `DB`;")
	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
# 	if (is.null(list$ref)==F){
# 		temp<-getcommon(temp)
# 	}
# 		print (temp[1:10,])
		table<-applynorm(temp[,c(1,2)],list)
		names(table)<-c('id',list$group1[ns])
# 		print(table[1:10,])
if (ns==1){
		all<-table
		names(all)<-c('id',list$group1[ns])
		
		#tableall<-temp
		
	}else{
		
		names(table)<-c('id',list$group1[ns])
		all<-merge(all,table,by="id",all=TRUE)
	
# 		templistnames<-
		#tableall<-merge(tableall,temp,all=TRUE,by="id")
	}
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1
# print (all[1:10,])
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
jpeg(paste(sep="",list$path,"img.jpg"))
plot(clust)
dev.off()
###############################################################
}else{
#####################################################
# print (table[1:10,])
d<-as.dist(1-cor((table[,2:ncol(table)]),method="spearman"))
# 
print(d)
clust<-hclust(d,method="complete")
# clust
jpeg(paste(sep="",list$path,"img.jpg"))
plclust(clust,labels=list$group1)
dev.off()
######################################################
}

doc = newXMLDoc()
top=newXMLNode("div",attrs=c(style="text-align: center;"))
newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = top)
save<-saveXML(top,file=paste(sep="",list$path,"result.html"))


dbDisconnect(con) 
