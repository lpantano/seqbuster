
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

if (is.null(list$freq1)==F){
	temp<-unlist(strsplit(list$freq1," "))
	options<-paste(sep=" ","where `freq` >",list$freq1," AND")
}else{
	options<-"where `freq` > 0 AND"
}
if (is.null(list$freq2)==F){
	temp<-unlist(strsplit(list$freq1," "))
	options<-paste(sep=" "," `freq` <",list$freq2," AND")
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

ns<-0
for (s in list$group1){
	names<-s
	#print(s)
	s<-paste(sep="",list$project,"`.`",s)
	
	ns<-ns+1
	
	query<-paste(sep="","select `id`,`freq` AS ROWS from `",s,"` ",options,";")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
# 	print (temp)
	if (ns==1){
		table<-temp
		names(table)<-c("id",names)
		
	}else{
		names(temp)<-c("id",names)
		table<-merge(table,temp,by="id",all=TRUE)
	}
}
otable<-table
table[is.na(table)]<-0
# table[1,]
if (list$norm=="log2"){
	table[,2:ncol(table)]<-log(table[,2:ncol(table)],base=2)
}else if (list$norm=="log10"){
	table[,2:ncol(table)]<-log(table[,2:ncol(table)],base=10)
}else if (list$norm=="ln"){
	table[,2:ncol(table)]<-log(table[,2:ncol(table)])
}
vector<-unlist(table[,2:ncol(table)])
#summary(vector[vector>0])
median<-quantile(vector,probs=.95)

scale1<-seq(0,median+1,by=(median+1)/5)
scale2<-seq(median+1,max(vector),by=(max(vector)-median)/5)
scale<-round(c(scale1,scale2[2:length(scale2)]))
# scale
values<-vector(length=9)
table2<-matrix(ncol=ncol(table)-1,nrow=9)
table3<-matrix(ncol=6,nrow=ncol(table)-1)

z<-1
for (i in 2:ncol(table)){
for (j in 1:9){
	table3[i-1,]<-summary(table[table[,i]>0,i])
	values[j]<-length(table[table[,i]>scale[j] & table[,i]<=scale[j+1] & table[,i]>0,i])
}
	#print(table[,j])
	table2[,z]<-values
	z<-z+1
	
}
otable<-table2
for (i in 1:ncol(table2)){
	table2[,i]<-table2[,i]/sum(table2[,i],na.rm=T)*100
}

# 
ttable<-t(table2)

col<-chosecolor(nrow(ttable))
# col<-col[1:nrow(ttable)]
namesscale<-paste(sep=":",format(scale[1:9],digits=3),format(scale[2:10],digits=3))
# namesscale
format(scale,digits=3,scientific=T)

jpeg(paste(sep="",list$path,"img.jpg"), quality=100)
 barplot(ttable,beside=T,col=col,names.arg=namesscale,legend.text=list$group1,xlab="Frequency range",ylab="Percentage of reads")
 #points(index[vector==T],maxvalue[vector==T]+2,pch="*",col="green",cex=2)
 dev.off()

namesrow<-c(namesscale)
namescol<-c("Fraction",list$group1)


library(XML)
doc = newXMLDoc()
top=newXMLNode("div",attrs=c(style="text-align: center;"))
newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = top)
tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)

tr<-newXMLNode("tr",  parent = tbl)
for (c in 1:(ncol(table2)+1)){
	td<-newXMLNode("td", attrs=c(bgcolor="yellow"),namescol[c], parent = tr)
# 	print (names(table[c]))
}

for (r in 1:nrow(table2))
{
tr<-newXMLNode("tr",  parent = tbl)

td<-newXMLNode("td", attrs=c(bgcolor="orange"),namesrow[r], parent = tr)
for (c in 1:ncol(table2))
{
	td<-newXMLNode("td", parent = tr)
	ac<-newXMLNode("acronym", attrs=c(title=otable[r,c]),round(table2[r,c],digit=2), parent = td)
	
}
}

newXMLNode("p",attrs=c(style= "font-weight:bold"),"Summary of the distribution", parent = top)

tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)

namestable3<-c("sample","min","Q25","median","mean","Q75","max")
tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
tr<-newXMLNode("tr",  parent = tbl)
for (c in 1:(ncol(table3)+1)){
	td<-newXMLNode("td", attrs=c(bgcolor="yellow"),namestable3[c], parent = tr)
# 	print (names(table[c]))
}
for (r in 1:nrow(table3))
{
tr<-newXMLNode("tr",  parent = tbl)

td<-newXMLNode("td", attrs=c(bgcolor="orange"),list$group1[r], parent = tr)
for (c in 1:ncol(table3))
{
	td<-newXMLNode("td", round(table3[r,c]),parent = tr)
# 	ac<-newXMLNode("acronym", attrs=c(title=otable[r,c]),round(table2[r,c],digit=2), parent = td)
	
}
}
saveXML(top,file=paste(sep="",list$path,"result.html"))
dbDisconnect(con) 

