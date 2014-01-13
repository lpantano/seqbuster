getgenes2<-function(data,algval,comalval,commival,pathval,projectval){

	################open micros
file<-read.table(data)
# file<-read.table(paste(sep="","/home/lpantano/micro/seqbuster/down"),header=TRUE)
# mirs<-"hsa-miR-34c-5p"
mifile<-as.vector(file[,1])
# print(mifile)
#########targets
common<-0
commonDB<-0
commonDB<-as.numeric(comalval)
common<-as.numeric(commival)
if (algval!=""){
	db<-c("targetscan","mirbase","pictar")
}else{
	db<-algval
	commonDB<-1
}
#print (list)
# tabletemp<-vector()
table<-data.frame(chr="0",nmir=0)
# nrow(file)
# query<-"select `GENEID` from `targetscan` where `miRNA` LIKE 'hsa-miR-379';"
# print(query)
# rs <- dbSendQuery(con,query) 
# temp<- as.data.frame(fetch(rs))
for (i in mifile){
# 	print(i)
	mirtarget<-data.frame(t=0,f=0)
	for (j in 1:length(db)){
		query<-paste(sep="","select `GENEID` from `",db[j],"` where `miRNA` LIKE '",i,"';")
		print(query)
		rs <- dbSendQuery(con,query) 
		temp<- as.data.frame(fetch(rs))
#  		print (nrow(temp))
		
		if (nrow(temp)>1){
			temp<-data.frame(t=unlist(unique(temp[,1])),f=rep(1,length(unlist(unique(temp[,1])))))
			names(temp)[2]<-db[j]
#  			print (nrow(temp))
			mirtarget<-merge(mirtarget,temp,all=T,by=1)
# 			print(mirtarget[1:10,])
			mirtarget[is.na(mirtarget)]<-0
		}
	}
# 	print(mirtarget[1:10,])
	if (length(mirtarget[mirtarget[,1]!=0,1])){
	mirtarget$sum<-rowSums(mirtarget[,2:ncol(mirtarget)])
# 	print(mirtarget[1:10,])
	mirtarget<-as.vector(unique(mirtarget[mirtarget$sum>=commonDB,1]))
	tabletemp<-data.frame(chr=mirtarget,nmir=rep(1,length(mirtarget)))
	table<-merge(table,tabletemp,by=1,all=T)
	table[is.na(table)]<-0
# 	print(table[1:10,])
	table$sum<-apply(table[,c(2,3)],1,sum)
# 	print(table[1:10,])
	table<-table[,c(1,4)]
	names(table)<-c("chr","nmir")
 }
# 	print(table[1:10,])
	
# 	target between micros	
}
table<-table[table$nmir>=common,]
 print(table[1:10,])
doc = newXMLDoc()
top=newXMLNode("html")
newXMLNode("p",parent=top)
# tbl<-newXMLNode("table",attrs=c(align="center"),parent=top)
# tr<-newXMLNode("tr",parent=tbl)
# col=lapply( names(targets),function(x) newXMLNode("th",x,parent=tr))
# addChildren(tr, col)
# if (nrow(targets)>0){
# for (i in 1:nrow(targets)){
# 	tr<-newXMLNode("tr",parsep = " ",ent=tbl)
# 	col=lapply(targets[i,],function(x) newXMLNode("td", x,parent=tr))
# 	addChildren(tr, col)
# }
# }
if (nrow(table)>0){
write.table(table,file=paste(sep="",pathval,"results.txt"),row.names = FALSE,col.names = FALSE,quote = FALSE,sep ="\t")
ref<-newXMLNode("a", attrs = c(href = paste(sep="","results.txt")),"Download list", parent = top)
#print ("Ok")
}else{
	p<-newXMLNode("p","no targets found with those parameters",parent=top)

}
saveXML(top,file=paste(sep="",pathval,"results.html"))

}
