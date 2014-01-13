list<-vector("list",length=3)
# list[[1]]<-c("mmlF16d","mml207d","mml1487d","mml3389d","mml9518d")
#   list[[1]]<-c("hsa2d","hsa204d","hsa5105d","hsa9277d","hsa19457d","hsa24090d","hsa32120d")
# list[[1]]<-c("dme01h","dme26h","dme610h")
 list[[1]]<-c("cel1L","cel2L","cel3L","cel4L")
# list[[1]]<-c("esc","eb")
namefile<-"cel.isos"
list[[3]]<-"t"
names(list)<-c("group1","g2","t")
list$project<-"spsmir"
list$DB<-"miRNA"
list$qmin<-0
list$qmax<-100
list$scale<-5000000
list$error<-"na"
list$trimmed5<-"trimmed5"
list$ref<-"ref"
list$path<-"/projects/srna_sps/data/"
library(RMySQL)
library(XML)
MySQL(max.con = 1, fetch.default.rec = 10000000, force.reload = FALSE)
source("/projects/NetBeansProjects/SeqBuster-2.0/Rscripts//R/stand.R")
m <- dbDriver("MySQL")
con <- dbConnect(m,host="localhost",user="lpantano",db="seqhand",password="sqllorena")


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
		
# 		mutation<-unlist(strsplit(gsub("[0-9]+","",m,""),""))
# 		nt1<-mutation[1]
# 		nt2<-mutation[2]
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

infogroup<-vector()
scale<-as.numeric(list$norm)
ns<-0
listsamples1<-vector("list",length=length(list$group1))
max1<-1:length(list$group1)
table<-data.frame(id=0,freq=0)
cof<-0

for (s in list$group1){
	#s<-paste(sep="",list$project,"`.`",s)
	s<-paste(sep="",s,list$project)
	ns<-ns+1
	#print(s)
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `id`,`seq`,`chr`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`DB`,`freq` AS freq from `",s,"` where `amb`=1 AND  `mut` like '0' ;")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))


	tempn<-applynorm(temp[,c(1,9)],list)
# 	print(tempn[1:10,])
	tempn<-merge(temp,tempn,by=1,all=FALSE)
	tempn<-tempn[tempn$DB==list$DB,]
	#tempn<-tempn[tempn$trimmed5!=0,]
#  	print(tempn[1:10,])
	temp<-tempn[,c(1,2,3,4,5,6,7,8,10)]
	names(temp)[9]<-"freq"
	temp$freq<-round(temp$freq)
#  	print(temp[1:10,])
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
# 				if (mi=="hsa-miR-124"){
# 					print (data.frame(id=listidtemp,freq=freqtemp));
# 				}
				table<-rbind(table,data.frame(id=listidtemp,freq=freqtemp))
# 				print (listidtemp)
		}else if(freqperfect>0 & is.null(list$ref)==F){
			labelper<-paste(sep="",mi,",Ref")
			table<-rbind(table,data.frame(id=labelper,freq=freqperfect))
		}

		
		
	}

	if (ns==1){
		all<-table[table$id!=0,]
		names(all)<-c('id',list$group1[ns])
		
	}else{
		
		names(table)<-c('id',list$group1[ns])
		all<-merge(all,table[table$id!=0,],by="id",all=TRUE)
	}
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1
table1<-all

table1<-table1[grep("tr5",table1$id),]
tableref<-all[grep("Ref",all$id),]

table<-table1[1:ns]
tableq<-table1[1:ns]

tabler<-tableref[1:ns]
tablerq<-tableref[1:ns]


for (i in 2:(ns)){

	table[,i]<-mapply(ztest,table1[,i],table1[,i+1],sum(table1[,i]),sum(table1[,i+1]))
	names(table)[i]<-paste(sep="-",names(table1)[i],names(table1)[i+1])
	sort<-sort(table[,i],index.return=T)
	ind<-1:nrow(table)
	order<-unlist(lapply(ind,function (x) ind[sort$ix==x]))
	tableq[,i]<-mapply(BHcorrection,table[,i],order,nrow(table))
	names(tableq)[i]<-paste(sep="-",names(table1)[i],names(table1)[i+1])
}

for (i in 2:(ns)){

	tabler[,i]<-mapply(ztest,tableref[,i],tableref[,i+1],sum(tableref[,i]),sum(tableref[,i+1]))
	names(tabler)[i]<-paste(sep="-",names(tableref)[i],names(tableref)[i+1])
	sort<-sort(tabler[,i],index.return=T)
	ind<-1:nrow(tabler)
	order<-unlist(lapply(ind,function (x) ind[sort$ix==x]))
	tablerq[,i]<-mapply(BHcorrection,tabler[,i],order,nrow(tabler))
	names(tablerq)[i]<-paste(sep="-",names(tableref)[i],names(tableref)[i+1])
}

#print only q$value <=0.05, some column
#row<- mir tr5 plot time serie, in each node text freq
meanrow<-rowSums(table1[2:ns+1])
indsort<-sort(meanrow,decreasing=T,index.return=T)

#scalecol<-c("#00CCFF","#0000FF","#006633","#99FF33","#CCFF00","#CC3333","#990066","BLACK")
scalecol<-(c("#FFCC99","#CC9966","#996633","#663300","#330000"))

#png(paste(sep="",list$path,"mml.iso.exp.png"),width = 520, height = 720, units = "px",res=200)
doc = newXMLDoc()
top=newXMLNode("body")
####par(mfrow=c(3,2))
count<-0
ten<-0
maxfreq<-0
minfreq<-1000000000
micros<-table1[1,]
tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
for (i in indsort$ix){
	count<-count+1
	#i<-indsort$ix[count]
	name<-unlist(strsplit(table[i,1],","))
	ref<-tableref[tableref$id==paste(sep="",name[1],",Ref"),]
	iso<-table1[i,]

	isov<-tableq[i,2:(ns)]
	refv<-tablerq[tableref$id==paste(sep="",name[1],",Ref"),2:(ns)]
	if(nrow(ref)==0){
			dif<-rep(1,ns)
			refv<-rep(1,ns-1)
			ref<-iso
			ref[1]<-paste(name[1],",Ref",sep="")
			ref[2:ns]<-refv
			ratio<-rep(100,ns)
	}else{
		dif<-refv[isov<0.05]
		total<-iso[2:ns]+ref[2:ns]
		ratio<-iso[2:ns]/total*100
	}
# 	if (length(isov[isov<0.05])){
# 		micros<-append(micros,table1[i,1])
# 	}
	if (length(dif)>0){
	if(max(dif)>0.05 & min(isov)<0.05){
	#ten<-ten+1
	colp<-as.character(cut(as.numeric(ratio),breaks=c(-1,20,40,60,80,101),labels=scalecol))
	tr<-newXMLNode("tr",  parent = tbl)
	td<-newXMLNode("td",table1[i,1] ,parent = tr)
	if (ten<5){
			

####		plot(1:ns,(table1[i,2:(ns+1)]),col=colp,bg=colp,pch=21,cex=1.5,ylab="expression RP5M",xlab="time series",main=table1[i,1],xaxt='n')
		if (maxfreq<max(log(table1[i,2:(ns+1)]))){
		maxfreq<-max(log(table1[i,2:(ns+1)]))
		}
		if (minfreq>max(log(table1[i,2:(ns+1)]))){
			minfreq<-min(log(table1[i,2:(ns+1)]))
		}
	}
	
	if (ten<5){
#####		axis(1,c(1:ns),labels=names(table1)[2:(ns+1)])
	}
	#coll<-as.character(cut(as.numeric(tableq[i,2:(ns)]),breaks=c(-1,0.051,2),labels=c("red","black")))
		one<-0
		td<-newXMLNode("td",table1[i,2] ,parent = tr)
		td<-newXMLNode("td",attrs=c(bgcolor="grey"),ref[2] ,parent = tr)
		for (j in 2:ns){
			coltr<-"white"
			coll<-"black"
			if (tableq[i,j]<0.05 & refv[j-1]>0.05){
				coll<-"red"
				coltr<-"red"
				one<-1
			}
			if (tableq[i,j]<0.05 & refv[j-1]<0.05){
				coll<-"grey"
				coltr<-"red"
			}
			td<-newXMLNode("td",attrs=c(bgcolor=coltr),table1[i,j+1] ,parent = tr)
			td<-newXMLNode("td",attrs=c(bgcolor="grey"),ref[j+1] ,parent = tr)
			if (ten<5){
						
####				lines((j-1):j,(table1[i,j:(j+1)]),col=coll)
			}
		}
		if (one==1) {
 			micros<-rbind(micros,iso)
			micros<-rbind(micros,ref)
		}
	}
	}
}

freq<-as.vector(as.matrix(all[,2:ns]))
maxtotal<-max(as.vector(as.matrix(all[,2:ns])))
colh<-c(rep("white",round(minfreq)),rep("blue",(round(maxfreq)-round(minfreq))),rep("white",(round(maxtotal)-round(minfreq))) )
#hist(log(freq[freq>1]),col=colh,main="frequency distribution of miRNAs/isomiRs",xlab="log units")
saveXML(top,file=paste(sep="",list$path,namefile,".html"))
micros$sum<-apply(micros[,2:ncol(micros)],1,sum)
write.table(micros[2:nrow(micros),], file=paste(sep="",list$path,namefile),row.names=F,quote=F,col.names=F,sep="\t")
dbDisconnect(con)
#dev.off()

###########################
# meanrow<-rowSums(tableref[2:ns])
# indsort<-sort(meanrow,decreasing=T,index.return=T)
# 
# for (i in indsort$ix[1:2]){
# 	minv<-min(tableref[i,2:(ns+1)])
# 	if (minv>0){
# 		##print (tableref[i,2:(ns+1)])
# 		name<-unlist(strsplit(tableref[i,1],","))
# 		print(name[1])
# 		iso<-table1[grep(paste(sep="",name[1],", tr5"),table1$id),]
# 		print(iso)
# 		if(nrow(iso)>0){
# 			minviso<-min(iso[2:(ns+1)])
# 			if (minviso==1){
# 				print (tableref[i,])
# 				print(iso)
# 			}
# 		}
# 	}
# 
# }
