
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
		if (pos>=as.numeric(list$ini) & pos<=as.numeric(list$size)){
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
if (is.null(list$DB)==F){
	
	options<-paste(sep="",options," `DB` LIKE '",list$DB,"' AND")
}
if (is.null(list$locus)==F){
	
	options<-paste(sep="",options," `chr` LIKE '%",list$locus,"%' AND")
}


temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])

cof<-as.numeric(list$cof)/100
# cof*2


# list$loci<-loci
jpeg(paste(sep="",list$path,"img.jpg"))
col<-rainbow(length(list$group1))
col=chosecolor(length(list$group1))
col<-col[1:length(list$group1)]
plot(1:10,1:10,type='n',ylim=c(0,25),xlim=c(0,2000),ylab="log(freq)",xlab="loci")
# plot(1,type='n')
legend("topright",legend=list$group1,fill=col)
scale<-as.numeric(list$norm)
ns<-0
# listsamples1<-vector("list",length=length(list$group1))
max1<-1:length(list$group1)
table<-data.frame(id=0,freq=0)
ma<-matrix(ncol=5,nrow=length(list$group1))
ma<-as.data.frame(ma)
data<-vector("list",length=length(list$group1))
#Load samples inside groups
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	#print(s)
	query<-paste(sep="","select `id`,`freq` from `",s,"` ORDER BY `freq` DESC;")
	print(query)
	rs <- dbSendQuery(con,query) 
	table <- as.data.frame(fetch(rs))
	max<-round(max(table[,2]),digits=3)
	total<-round(sum(table[,2]),digits=3)
	n<-nrow(table)
	r<-round(max/total,digits=3)
	ma[ns,]<-c(s,n,total,max,r)
	table<-applynorm(table,list)
	points(sort((table[,2]),decreasing=TRUE),col=col[ns],cex=.2)
	data[[ns]]<-unlist(unique(table[,2]))
}
dev.off()

cor<-matrix(ncol=length(list$group1),nrow=length(list$group1))
for (i in 1:length(list$group1)){
for (j in 1:length(list$group1)){

	wc<-wilcox.test(data[[i]],data[[j]],paired=FALSE)
	cor[i,j]<-wc$p.value

}
}

# print (cor)


doc = newXMLDoc()
top=newXMLNode("div",attrs=c(style="text-align: center;"))
newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = top)
table<-ma
names(table)<-c("sample","reads","total","max","ratio")
tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)

tr<-newXMLNode("tr",  parent = tbl)

for (c in 1:ncol(table)){
	td<-newXMLNode("td", attrs=c(bgcolor="yellow"),names(table[c]), parent = tr)
# 	print (names(table[c]))
}

for (r in 1:nrow(table)){
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(table[r,],function(x) newXMLNode("td", x,parent=tr))
	addChildren(tr, col)
	
}
p<-newXMLNode("p",parent=top)
#correlation
p<-newXMLNode("p","similarity between samples is calculated with
the wilcoxon test giving a p.value for each pair of samples. Values
with a p.value lower than 0.05 indicate differente experiment
capacities between samples.",parent=top)
tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)

tr<-newXMLNode("tr",  parent = tbl)
td<-newXMLNode("td","Samples", parent = tr)
for (c in 1:ncol(cor)){
	td<-newXMLNode("td", attrs=c(bgcolor="orange"),list$group1[c], parent = tr)
# 	print (names(table[c]))
}

for (r in 1:nrow(cor)){
	tr<-newXMLNode("tr",  parent = tbl)
	td<-newXMLNode("td", attrs=c(bgcolor="yellow"),list$group1[r], parent = tr)
	col=lapply(round(cor[r,],digits=3),function(x) newXMLNode("td", x,parent=tr))
	addChildren(tr, col)
	
}

save<-saveXML(top,file=paste(sep="",list$path,"result.html"))
dbDisconnect(con) 
