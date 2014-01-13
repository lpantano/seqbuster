
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

parsename<-function(name){
	temp<-unlist(strsplit(name,","))
	temp<-unlist(strsplit(temp[2],":"))
	if (length(temp)==1){
		temp[1]<-"NA"
	}
	if (length(grep("Ref",name))>0){
		temp[1]<-"Ref"
	}
#  	print(temp[1])
	return(temp[1])
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
if (is.null(list$DB)==F){
	
	options<-paste(sep="",options," `DB` LIKE '",list$DB,"' AND")
}
if (is.null(list$locus)==F){
	
	options<-paste(sep="",options," `chr` LIKE '%",list$locus,"%' AND")
}
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
if (is.null(list$error)){
	list$error<-"0"

}
#Load samples inside groups
for (s in list$group1){
	sname<-s
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	#print(s)
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `id`,`seq`,`chr`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`DB`,`freq` AS freq from `",s,"` ",options," AND `chr` NOT LIKE 'amb' and `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+' AND `mut` NOT REGEXP '^[0-9]+[ACTG]+[0-9]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	
# 	if (is.null(list$ref)==F){
# 		temp<-getcommon(temp)
# 	}
	
	tempn<-applynorm(temp[,c(1,9)],list)
# 	print(tempn[1:10,])
	tempn<-merge(temp,tempn,by=1,all=FALSE)
# 	print(tempn[1:10,])
	temp<-tempn[tempn$DB==list$DB,c(1,2,3,4,5,6,7,8,10)]
	names(temp)[9]<-"freq"
	temp$freq<-round(temp$freq)
# 	print(temp[1:10,])
	query<-paste(sep="","select SUM(freq) AS freq from `",s,"`;")
	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
	table<-data.frame(id=0,freq=0)
	listmitemp<-unlist(unique(temp$chr))
	for (mi in listmitemp){
	
		pertemp<-temp[temp$chr==mi & temp$trimmed5==0 & temp$trimmed3==0 & temp$mut==0 & temp$addition3==0,]
		freqperfect<-0
		if (length(pertemp$chr)>0){freqperfect<-pertemp$freq}
		vartemp<-temp[temp$chr==mi & (temp$trimmed5!=0 | temp$trimmed3!=0 | temp$mut!=0 | temp$addition3!=0),]
		tempvalue<-sum(vartemp$freq)
		
		minvalue<-min(freqperfect,tempvalue,na.rm=T)
		
		vartemp<-filter(vartemp,minvalue*cof,7,type="filter",freq=freqperfect,error=list$error)
		
		
		if (length(vartemp$id)>0){
# 			for (nvar in 1:length(vartemp$id)){
# 				totalfreqid<-sum(vartemp$freq)
				#get all freq variation and cutoff 10%	
# 				cutofftemp<-cof*totalfreqid
# 				vartemp<-vartemp[vartemp$freq>=cutofftemp,]
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
# 		if (mi=="hsa-miR-22"){
# 			print (table)
# # 			print (listidtemp)
# # 			print (freqtemp)
# # 			print("kk")
# 		}
		
		
	}
# 	print (nrow(table))
#	listsamples1[[ns]]<-temp
# 	table$n<-round(table[,2])
	if (ns==1){
		all<-table[table$id!=0,]
		names(all)<-c('id',sname)
		#tableall<-temp
		
	}else{
		
		names(table)<-c('id',sname)
		all<-merge(all,table[table$id!=0,],by="id",all=TRUE)
		#tableall<-merge(tableall,temp,all=TRUE,by="id")
	}
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1
table2<-all
# print(table2[1:10,])
table2$type<-as.character(mapply(parsename,table2[,1]))
#print(table2)
# table2
types<-unlist(unique(table2$type))
col<-chosecolor(length(list$group1))
# table2<-table2[order(table2[,2],decreasing=TRUE),]
jpeg(paste(sep="",list$path,"img.jpg"),width=500,height=(length(types))*450,units='px')
par(mfrow=c(length(types),1),mar=c(4,4,4,4))

nnt<-0
error<-0
for (nt in types){
	nnt<-nnt+1
	#print(nt)
	tabletemp<-table2[table2$type==as.character(nt),]
#  	print(tabletemp)
	if (nrow(tabletemp)<10 | length(list$group1)!=2){
		error<-"not enough observations or you have only selected one sample."

	}else{
# 	print(nrow(tabletemp))
	corvalue<-cor.test((tabletemp[,2]),(tabletemp[,3]))
	# corvalue
	coltemp<-col[nnt]
	plot(log(tabletemp[,2]),log(tabletemp[,3]),type='p',ylab=list$group1[2],xlab=list$group1[1],cex=0.2,sub=paste("r: ",round(corvalue$estimate,digits=2)),main=paste("type of sequence: ",nt))
	}
}
dev.off()

doc = newXMLDoc()
top=newXMLNode("div",attrs=c(style="text-align: center;"))
newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = top)

if (error!=0){
	newXMLNode("b",error,parent=top)
}
saveXML(top,file=paste(sep="",list$path,"result.html"))
dbDisconnect(con) 
