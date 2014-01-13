nochanges<-function(len1val,len2val,freq1val,freq2val,samplesnames,pathval,dbval,projectval){
	if (freq1val!=""){
		options<-paste(sep=" ","where `freq` ",freq1val," AND ")
	}else{
		options<-"where `freq` > 0 AND "
	}
	if (freq2val!=""){
		options<-paste(sep=" ",options," `freq` ",freq2val," AND ")
	}
	
	if (len1val!=""){
		
		options<-paste(sep=" ",options,"`len` ",len1val," AND")
	}
	if (len2val!=""){
		
		options<-paste(sep=" ",options,"`len` ",len2val," AND")
	}
	if (dbval!=""){
		
		options<-paste(sep="",options," `DB` LIKE '",dbval,"' AND")
	}
	
temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])


nts<-c("A","T","C","G")
ntspos<-c(1,2,3,4)

ns<-1
listsamples<-vector("list",length=length(samplesnames))
lensamples<-0
freq<-vector()
names<-samplesnames
for (s in samplesnames){
	#s<-paste(sep="",s,projectval)
	query<-paste(sep="","select `chr`,`freq` from `",projectval,"`.`",s,"` ",options," AND `chr` NOT LIKE 'amb' AND  `trimmed5` like '0' AND  `trimmed3` like '0' AND  `addition3` like '0' AND  `mut` like '0' group by `name`;")
# 	print(query)
	rs <- dbSendQuery(con,query) 
	perf <- as.data.frame(fetch(rs))
# 	print(perf[,1])
	query<-paste(sep="","select `chr`,`freq` from `",projectval,"`.`",s,"` ",options," AND `chr` NOT LIKE 'amb' AND  (`trimmed5` not like '0' OR  `trimmed3` not like '0' OR  `addition3` not like '0' OR  `mut` not like '0') group by `name`;")
	
	rs <- dbSendQuery(con,query) 
	all <- as.data.frame(fetch(rs))
	
	tablemirs<-merge(perf,all,by=1,all=T)
	tablemirs[is.na(tablemirs)]<-0
# 	print(tablemirs[1:10,])
	temp<-tablemirs[tablemirs[,3]==0,c(1,2)]
# 	temp<-data.frame(name=nochange,type=rep(1,length(nochange)))
	names(temp)<-c("name",s)
	if (ns==1){
		listmi<-temp
	}else{
		listmi<-merge(temp,listmi,by=1,all=TRUE)
	}
	ns<-ns+1
	freq<-unlist(unique(c(freq,temp$freq)))
}

listmi[is.na(listmi)]<-0
# listmi[1:10,]

doc = newXMLDoc()
top=newXMLNode("html")
script<-newXMLNode("script", attrs=c(language="JavaScript",src=paste(sep="",pathjava,"javascript/showtableid.js")),"",parent=top)

tblbig<-newXMLNode("table",attrs=c(align="center",border=1,height="10px",width="20%"),parent=top)
trbig<-newXMLNode("tr",parent=tblbig)
col=lapply(c("locus",names),function(x) newXMLNode("th",attrs=c(bgcolor="orange"),x,parent=trbig))
addChildren(trbig, col)


for (i in 1:nrow(listmi)){
	trbig<-newXMLNode("tr",parent=tblbig)
	color<-"white"
	if (prod(listmi[i,2:ncol(listmi)])>0){
		color<-"#FF6666"
	}
	trbig<-newXMLNode("td",attrs=c(bgcolor=color),as.character(listmi[i,1]),parent=trbig)
	text<-unlist(listmi[i,2:ncol(listmi)])
# 	print(text)
	col=lapply(text,function(x) newXMLNode("td",x,parent=trbig))
	addChildren(trbig, col)

}

saveXML(top,file=paste(sep="",pathval,"results.html"))

}
