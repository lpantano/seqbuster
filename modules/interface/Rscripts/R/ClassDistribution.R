
getind<-function(n,all){
	vector<-vector()
	num<-1:length(all)
	for (i in n){
		if(length(all[all==i])>0){
			ind<-num[all==i]
			vector<-append(vector,ind)
		}
		
	}
	return(vector)
}


library("XML")
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")
lp<-parse("temp/script_param")


MySQL(max.con = 1, fetch.default.rec = 10000000, force.reload = FALSE)
m <- dbDriver("MySQL")
if (port==0){
	con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw)
}else{
	con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw,port=port)
}

if (is.null(lp$freq1)==F){
	temp<-unlist(strsplit(lp$freq1," "))
	options<-paste(sep=" ","where `freq` <=",lp$freq1," AND")
}else{
	options<-"where `freq` > 0 AND"
}
if (is.null(lp$freq2)==F){
	temp<-unlist(strsplit(lp$freq1," "))
	options<-paste(sep=" "," `freq` >=",lp$freq2," AND")
}
if (is.null(lp$len1)==F){
	temp<-unlist(strsplit(lp$len1," "))
	options<-paste(sep=" ","  `len` <=",lp$len1," AND")
}
if (is.null(lp$len2)==F){
	temp<-unlist(strsplit(lp$len2," "))
	options<-paste(sep=" "," `len` >=",lp$len2," AND")
}
if (is.null(lp$linker)==F){
	
	options<-paste(sep="",options," `tag` LIKE '",lp$linker,"' AND")
}
if (is.null(lp$DB)==F){
	
	options<-paste(sep="",options," `DB` LIKE '",lp$DB,"' AND")
}
if (is.null(lp$locus)==F){
	
	options<-paste(sep="",options," `chr` LIKE '%",lp$locus,"%' AND")
}

#options

temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])

counter<-"(freq"
if (lp$operation=="Different locus"){
	
	operation<-"COUNT(DISTINCT "
	counter<-"chr"
}else if (lp$operation=="Total frequency"){
	operation<-"sum"
}else if (lp$operation=="Different reads"){
	operation<-"COUNT"
}else if (lp$operation=="MAX"){
	operation<-"max"
}else if (lp$operation=="MIN"){
	operation<-"min"
}



if (lp$dis=="DB"){
	dist<-"DB"
}else{
	dist<-"len"
}

ns<-0
for (s in lp$group1){
	ns<-ns+1
	#print(s)
	sname<-s
	s<-paste(sep="",lp$project,"`.`",s)
	query<-paste(sep="","select `",dist,"`, ",operation,counter,") AS ROWS from `",s,"` ",options," group by `",dist,"`;")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
# 	print (temp)
	if (ns==1){
		table<-temp
		names(table)<-c(dist,sname)
		
	}else{
		names(temp)<-c(dist,sname)
		table<-merge(table,temp,by=dist,all=TRUE)
	}
}
otable<-table

cols<-length(lp$group1)+1
table[is.na(table)]<-0
#normilize
if (lp$norm=="1"){
for (i in 2:cols){
	table[,i]<-table[,i]/sum(table[,i],na.rm=T)*100
}
}else{
	for (i in 2:cols){
		table[,i]<-table[,i]
	}

}
# table

if (lp$tchart=="Bars chart"){
jpeg(paste(sep="",lp$path,"img.jpg"))

ttable<-t(table[,2:ncol(table)])
col<-chosecolor(ncol(table)-1)
# jpeg(paste(sep="",path,"pictures/",args[1],".bmp"))
barplot(ttable,beside=T,col=col,names.arg=table[,1],legend.text=lp$group1)
#points(index[vector==T],maxvalue[vector==T]+2,pch="*",col="green",cex=2)
}else{
	jpeg(paste(sep="",lp$path,"img.jpg"),width=10,height=length(lp$group1)*7)
	
	numchart<-length(lp$group1)+2
# 	if (length(samplesnames)/2==round(length(samplesnames)/2)){
# 		numchart<-numchart+1
# 	}
	par(mfrow=c(numchart,2),plt=c(0.15,0.85,0.05,0.85))
	namedb<-""
	col<-vector()
	for (i in 2:(ncol(table))){
		cutoff<-sum(table[,i])*0.025
		temp<-table[table[,i]>=cutoff,c(1,i)]
		
		temp<-rbind(temp,c("Others",sum(table[table[,i]<cutoff,i])))
# 		print(temp)
		namesdbtemp<-table[table[,i]>=cutoff,1]
# 		print(namesdbtemp)
# 		print(namedb)
		difnames<-setdiff(namesdbtemp,namedb)
		
# 		print(difnames)
		if (length(difnames)>0){
			namedb<-c(namedb,difnames)
			coltemp<-chosecolor(length(namedb))
			col<-append(col,coltemp[(length(namedb)-length(difnames)+1):length(namedb)])
		}
		
# 		n<-nrow(temp)
# 		print(col)
		indcol<-getind(namesdbtemp,namedb)
		colvector<-col[indcol]
# 		col<-chosecolor(n-1)
		pie(as.numeric(temp[,2]),col=c(colvector,"grey"),labels="",main=names(table[i]))
		plot(1,type="n",axes=F,xlab="", ylab="")
		legend("topleft",c(table[table[,i]>=cutoff,1],"Others"),fill=c(colvector,"grey"))
		
	}
	

}
dev.off()


doc = newXMLDoc()
top=newXMLNode("div",attrs=c(style="text-align: center;"))
newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = top)
tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)

tr<-newXMLNode("tr",  parent = tbl)
for (c in 1:ncol(table)){
	td<-newXMLNode("td", attrs=c(bgcolor="yellow"),names(table[c]), parent = tr)
	print (names(table[c]))
}

for (r in 1:nrow(table))
{
tr<-newXMLNode("tr",  parent = tbl)

td<-newXMLNode("td", attrs=c(bgcolor="orange"),table[r,1], parent = tr)
for (c in 2:ncol(table))
{
	td<-newXMLNode("td", parent = tr)
	ac<-newXMLNode("acronym", attrs=c(title=otable[r,c]),round(table[r,c],digit=2), parent = td)
	
}
}

dbDisconnect(con)
saveXML(top,file=paste(sep="",lp$path,"result.html"))






