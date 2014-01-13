
library("XML")
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")
list<-parse("temp/script_param")


MySQL(max.con = 1, fetch.default.rec = 10000000, force.reload = FALSE)
m <- dbDriver("MySQL")
con <- dbConnect(m,host=hostname,user=username,db=dbname,password=pssw)


splitadapter<-function (str){

	
# 	print(vector)
	adapsummary<-matrix(ncol=10,nrow=4)
	adapsummary[is.na(adapsummary)]<-0
# 	adap<-matrix(ncol=40,nrow=nrow(temp))
	for (nr in 1:nrow(str)){
# 	print(adap)
	
		vector<-as.numeric(unlist(strsplit(str[nr,2]," ")))
# 		ind<-seq(1,length(vector),2)
# 		print (vector)
# 		vector<-unlist(vector[ind])
		
		for (i in 1:10){
	# 		print(vector)
			adapsummary[1,i]<-(vector[i])+adapsummary[1,i]
			adapsummary[2,i]<-(vector[i+10*1])+adapsummary[2,i]
			adapsummary[3,i]<-(vector[i+10*2])+adapsummary[3,i]
			adapsummary[4,i]<-(vector[i+10*3])+adapsummary[4,i]
		}
	# 	return((vector))
	}
	adapsummary<-scale(adapsummary,center=FALSE,scale=colSums(adapsummary))
	adapsummary[is.na(adapsummary)]<-0
	return(adapsummary)
}
# neededop<-list(mut="",trimmed3="",trimmed5="",addition3="")

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
if (is.null(list$seq)==F){
	
	options<-paste(sep="",options," `seq` LIKE '%",list$seq,"%' AND")
}

temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])

listsamples<-vector("list",length=length(list$group1))
adap<-"infoadap5"
if (list$adapter=="3 adapter"){
	adap<-"infoadap3"
}
ns<-0
lensamples<-length(list$group1)
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	ns<-ns+1
	query<-paste(sep="","select `id`,`",adap,"` from `",s,"` ",options,";")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
# 	print(temp)
# 	len<-length(splitadapter(temp[1,2]))
	listsamples[[ns]]<-splitadapter(temp)
#  	listsamples[[ns]]<-temp
	
}

listmaN<-listsamples
# listsamples[[1]]
listypos<-vector("list",(lensamples))
listxpos<-vector("list",(lensamples))
summaryN<-matrix(ncol=10,nrow=(lensamples))

# 	maxsummary<-matrix(ncol=ncol(listmaN[[1]]),nrow=length(lensamples))
for (i in 1:lensamples){

	listmaN[[i]]<-listmaN[[i]]*100
	listxpos[[i]]<-seq(i,i+((lensamples)+2)*(ncol(listmaN[[1]]))-1,(
lensamples)+2)
	listypos[[i]]<-listmaN[[i]]
	listypos[[i]][1,]<-listmaN[[i]][1,]
	listypos[[i]][2,]<-apply(listmaN[[i]][1:2,],2,sum)
	listypos[[i]][3,]<-apply(listmaN[[i]][1:3,],2,sum)
	listypos[[i]][4,]<-apply(listmaN[[i]][1:4,],2,sum)
		
	summaryN[i,]<-apply(listmaN[[i]][1:4,],2,sum)
		
}	

maxyplot<-max(summaryN)+max(summaryN)*0.20
maxxplot<-((lensamples)+1)*14+2
coln<-c("red","blue","green","yellow")
widthpic<-10*3
# 	bitmap(paste(sep="","/srv/www/htdocs/demo2/pictures/",type,".bmp"),height = 30, width = widthpic,units='cm',type="png16m")
bitmap(paste(sep="",list$path,"img.bmp"),height = 20, width =widthpic ,units='cm',type="png16m")
plot(1, type="n", axes=F, xlab="positions", ylab="",xlim=c(-1,maxxplot+1),ylim=c(0,maxyplot+0.20*maxyplot),main="Adapter quality",cex=2)
legend("top",c("A","U","C","G"),fill=coln, cex=1,horiz=TRUE)
#xaxis
segments(0,0,maxxplot,0)
#yaxis
segments(0,maxyplot,0,0)
yaxes<-seq(10,maxyplot,10)
segments(0,yaxes,-0.2,yaxes)
text(maxxplot,maxyplot/2,"% of NTs",srt=90,cex=.7)
text(-0.1,yaxes,labels=yaxes,pos=2,cex=.7)
#xlabs
text((listxpos[[1]]+listxpos[[(lensamples)]]+1)/2,0,1:10,cex=1,pos=1)
#print bars
for (i in 1:lensamples){
# 		print(listma[[i]])
		#micros
# 		pvect<-0
	rect(listxpos[[i]],0,listxpos[[i]]+1,listypos[[i]][1,],col=coln[1])
	rect(listxpos[[i]],listypos[[i]][1,],listxpos[[i]]+1,listypos[[i]][2,],col=coln[2])
	rect(listxpos[[i]],listypos[[i]][2,],listxpos[[i]]+1,listypos[[i]][3,],col=coln[3])
	rect(listxpos[[i]],listypos[[i]][3,],listxpos[[i]]+1,listypos[[i]][4,],col=coln[4])
		
}
dev.off()

# print(summary)
doc = newXMLDoc()
top=newXMLNode("div",attrs=c(style="text-align: center;"))
newXMLNode("img", attrs = c(src = paste(sep="","img.bmp")), parent = top)
# newXMLNode("p",paste("Correlation: ",corvalue$estimate), parent = top)
saveXML(top,file=paste(sep="",list$path,"result.html"))

dbDisconnect(con) 
