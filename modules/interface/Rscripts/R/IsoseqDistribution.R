library("XML")
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")
list<-parse("temp/script_param")

getvalueadd<-function(now,add,total,opt){
	value<-now
	if (opt=="ALL"){
		value<-c(now,add)
	}else if(opt=="permut"){
		if(add>=(0.1*total)){
			value<-c(now,add)
# 			print ("ok")
		}
	}
	return (value)
}

getvalue<-function(now,opt){
	if (opt=="MEDIAN"){
		value<-median(now)
	}else if(opt=="MEAN"){
		value<-mean(now)
	}else if(opt=="3Q"){
		value<-quantile(now,probs=.75)
	}else if(opt=="MAX"){
		value<-max(now)
	}else if(opt=="num"){
		value<-length(now)
	}
	return (value)
}



writexml<-function(listma,lensamples,names,labelsmut){
	listmaN<-listma
	#print (listma)
	listxpos<-vector("list",length(lensamples))
	listypos<-vector("list",length(lensamples))
	summary<-matrix(ncol=ncol(listmaN[[1]]),nrow=length(lensamples))
	summaryN<-matrix(ncol=ncol(listmaN[[1]]),nrow=length(lensamples))
# 	maxsummary<-matrix(ncol=ncol(listmaN[[1]]),nrow=length(lensamples))
	for (i in 1:length(lensamples)){
		listmaN[[i]]<-listma[[i]]/lensamples[i]*100
		listxpos[[i]]<-seq(i,i+(length(lensamples)+2)*(ncol(listmaN[[1]]))-1,length(lensamples)+2)
		listypos[[i]]<-listmaN[[i]]
		if (ncol(listmaN[[i]])==1){
			listypos[[i]][1,]<-listmaN[[i]][1,1]
			listypos[[i]][2,]<-sum(listmaN[[i]][1:2,1])
			listypos[[i]][3,]<-sum(listmaN[[i]][1:3,1])
			listypos[[i]][4,]<--(listmaN[[i]][4,1])
			listypos[[i]][5,]<--sum(listmaN[[i]][4:5,1])
			listypos[[i]][6,]<--sum(listmaN[[i]][4:6,1])
			listypos[[i]][7,]<--sum(listmaN[[i]][4:7,1])
			listypos[[i]][8,]<--sum(listmaN[[i]][4:8,1])
			summary[i,]<-sum(listma[[i]][1:3,1])
			summaryN[i,]<-sum(listmaN[[i]][1:3,1])
		
		
		}else{
			listypos[[i]][1,]<-listmaN[[i]][1,]
			listypos[[i]][2,]<-colSums(listmaN[[i]][1:2,])
			listypos[[i]][3,]<-colSums(listmaN[[i]][1:3,])
			listypos[[i]][4,]<--(listmaN[[i]][4,])
			listypos[[i]][5,]<--colSums(listmaN[[i]][4:5,])
			listypos[[i]][6,]<--colSums(listmaN[[i]][4:6,])
			listypos[[i]][7,]<--colSums(listmaN[[i]][4:7,])
			listypos[[i]][8,]<--colSums(listmaN[[i]][4:8,])
			summary[i,]<-colSums(listma[[i]][1:3,])
			summaryN[i,]<-colSums(listmaN[[i]][1:3,])
		}
	}	
	
	maxyplot<-100
	maxxplot<-(length(lensamples)+2)*ncol(listma[[1]])
	coln<-c("white","grey","black")
# 	colr<-(c("#151B8D","#6698FF","#347C2C","#CCFB5D","#FFFF00"))
	colr<-rev(c("#FFCC99","#CC9966","#996633","#663300","#330000"))
	colall<-c(coln,colr)
	jpeg(paste(sep="",list$path,"img.jpg"),quality=100)
	par(mar=c(1,1,1,1))
	plot(1, type="n", axes=F, xlab="", ylab="",xlim=c(-1,maxxplot+1),ylim=c(-maxyplot-0.40*maxyplot,maxyplot+0.20*maxyplot),main="Type of variability")
	legend("top",c("1 change","2 changes","+2 changes"),fill=coln, cex=.7,horiz=TRUE)
	#xaxis
	segments(0,0,maxxplot,0)
	#yaxis
	segments(0,maxyplot,0,0)
	yaxes<-seq(10,maxyplot,10)
	segments(0,yaxes,-0.2,yaxes)
	text(maxxplot,maxyplot/2,"% of miRNAs",srt=90,cex=.7)
	text(-0.2,yaxes,labels=yaxes,pos=2,cex=.7)
	#yaxisminus
	segments(maxxplot,-maxyplot,maxxplot,0)
	yaxes<-seq(-10,-maxyplot,-10)
	segments(maxxplot,yaxes,maxxplot+0.2,yaxes)
	text(-0.2,-maxyplot/2,"% of miRNAs",srt=90,cex=.7)
	text(maxxplot+0.2,yaxes,labels=-yaxes,pos=4,cex=.7)
	#legend
	poslegend<-seq(maxxplot/2-0.2*maxxplot,length(colr)*0.05*maxxplot+maxxplot/2-0.2*maxxplot,0.05*maxxplot)
	yvar<-0.07*maxyplot
	ylegend<-seq(-maxyplot-0.20*maxyplot,-maxyplot-0.20*maxyplot+yvar,1)
# 	
	polygon(c(poslegend[1],poslegend[2],poslegend[2]),c(ylegend[1],ylegend[1],ylegend[2]),col=colr[5],border='NA')
	polygon(c(poslegend[2],poslegend[2],poslegend[3],poslegend[3]),c(ylegend[1],ylegend[2],ylegend[3],ylegend[1]),col=colr[4],border='NA')
	polygon(c(poslegend[3],poslegend[3],poslegend[4],poslegend[4]),c(ylegend[1],ylegend[3],ylegend[4],ylegend[1]),col=colr[3],border='NA')
	polygon(c(poslegend[4],poslegend[4],poslegend[5],poslegend[5]),c(ylegend[1],ylegend[4],ylegend[5],ylegend[1]),col=colr[2],border='NA')
	polygon(c(poslegend[5],poslegend[5],poslegend[6],poslegend[6]),c(ylegend[1],ylegend[5],ylegend[6],ylegend[1]),col=colr[1],border='NA')
	text(mean(poslegend),-0.20*maxyplot-maxyplot,"variability significance",pos=1,cex=0.7)
# 	rect(poslegend,ylegend,poslegend+2,ylegend-0.05*maxyplot,col=colr)
	#group of mutation
	ylegend<--maxyplot
	text((listxpos[[1]]+listxpos[[length(lensamples)]]+1)/2,ylegend,labelsmut[2:length(labelsmut)],cex=1,pos=1)
	
# 	rect(poslegend,ylegend,poslegend+1,ylegend-1,col=col)
# 	print(listma)
	namesind<-1:length(names)
	for (i in 1:length(lensamples)){
# 		print(listxpos[[i]])
		#####samples labels
		text(listxpos[[i]]+0.5,listypos[[i]][3,],namesind[i],cex=.7,pos=3)
		#micros
		rect(listxpos[[i]],0,listxpos[[i]]+1,listypos[[i]][1,],col=coln[1])
		rect(listxpos[[i]],listmaN[[i]][1,],listxpos[[i]]+1,listypos[[i]][2,],col=coln[2])
		rect(listxpos[[i]],listypos[[i]][2,],listxpos[[i]]+1,listypos[[i]][3,],col=coln[3])
		#ratios
		rect(listxpos[[i]],0,listxpos[[i]]+1,listypos[[i]][4,],col=colr[1])
		rect(listxpos[[i]],listypos[[i]][4,],listxpos[[i]]+1,listypos[[i]][5,],col=colr[2])
		rect(listxpos[[i]],listypos[[i]][5,],listxpos[[i]]+1,listypos[[i]][6,],col=colr[3])
		rect(listxpos[[i]],listypos[[i]][6,],listxpos[[i]]+1,listypos[[i]][7,],col=colr[4])
		rect(listxpos[[i]],listypos[[i]][7,],listxpos[[i]]+1,listypos[[i]][8,],col=colr[5])
	}
	ytop<-maxyplot
	var<-0.05*maxyplot
	colp<-c("red","blue","green","yellow")
	nc<-1
	if (is.null(list$pvaluetype)==F){
		for (m1 in 1:(ncol(summary) ) ) {
			
			nc<-nc+1
			for (m2 in (m1):ncol(summary)){
				p<-pvalue(mean(lensamples),mean(lensamples),mean(summary[,m1]),mean(summary[,m2]))
				if(p<0.05){
					ytop<-ytop+var
					x1<-(listxpos[[1]][m1]+listxpos[[length(lensamples)]][m1]+1)/2
					x2<-(listxpos[[1]][m2]+listxpos[[length(lensamples)]][m2]+1)/2
					y1<-max(summaryN[,m1])+var
					y2<-max(summaryN[,m2])+var
	# 				ytop<-max(summary)
					segments(x1,y1,x1,ytop,colp[nc])
					segments(listxpos[[1]][m1],y1,listxpos[[length(lensamples)]][m1]+1,y1,colp[nc])
					segments(x2,y2,x2,ytop,colp[nc])
					segments(listxpos[[1]][m2],y2,listxpos[[length(lensamples)]][m2]+1,y2,colp[nc])
					segments(x1,ytop,x2,ytop,colp[nc])
	# 				print(c(m1,m2,p))
				}
			}
		}
	}
	dev.off()
	
	doc = newXMLDoc()
	top=newXMLNode("html")
	div<-newXMLNode("div",attrs=c(style="text-align: center;"),parent=top)
	newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = div)
	script<-newXMLNode("script", attrs=c(language="JavaScript",src="showtableid.js"),"",parent=top)
	tblbig<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
	trbig<-newXMLNode("tr",  parent = tblbig)
	tdbig<-newXMLNode("td","lable" ,parent = trbig)
	tdbig<-newXMLNode("td","name" ,parent = trbig)
	for (ind in namesind){
		trbig<-newXMLNode("tr",  parent = tblbig)
		tdbig<-newXMLNode("td",ind ,parent = trbig)
		tdbig<-newXMLNode("td",names[ind] ,parent = trbig)
	}
	

# 	colrnames<-c("0-20","20-40","40-60","60-80","80-100")
	colrnames<-c(">80","60-80","60-40","40-20","<20")
	br<-newXMLNode("p","Percentage of the isomiRs with respect to the corresponding reference miRNAs. ",attrs=c(align="center"),parent=top)
	
	text<-"The variability significance is calculated using the following equation: r=Fv/(Fr+Fv)*100, where Fr is the frequency of the reference sequence and Fv is the frequency of the variant sequence."
	newXMLNode("a","?",attrs=c(title=text,href=""),parent=br)
	tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply((colr),function(x) newXMLNode("td", attrs=c(bgcolor=x,style="height:1em;"),parent=tr))
	addChildren(tr, col)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(colrnames,function(x) newXMLNode("td", x,parent=tr))
	addChildren(tr, col)

	newXMLNode("p",parent=top)
		
	labelsrow<-c("ONE VARIANT","TWO VARIANTS","MORE","+80","80-60","60-40","40-20","20-0")
	tblbig<-newXMLNode("table",attrs=c(align= "center",border=0), parent = top)
	trbig<-newXMLNode("tr",  parent = tblbig)
	for (nm in 2:length(labelsmut)){
		tdbig<-newXMLNode("td",  parent = trbig)
		tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = tdbig)
		tr<-newXMLNode("tr",  parent = tbl)
			
# 		for (ns in 1:length(names)){
			td<-newXMLNode("td",attrs=c(bgcolor="orange"),labelsmut[nm], parent = tr)
# 		}
		for (nr in 1:8){
			tr<-newXMLNode("tr",  parent = tbl)
			
# 			for (ns in 1:length(names)){
				index<-ns*8-8+nr
				
# 				text<-paste(listmu[[nm]][[index]])
				#text<-paste(sep="","#",nm,nr)	
				text<-paste(sep="","javascript:loadtable('",nm,nr,"')")			
				td<-newXMLNode("td",attrs=c(bgcolor=colall[nr]),parent = tr)
				a<-newXMLNode("a",attrs=c(href=text),labelsrow[nr],parent=td)
# 			}
		}
	
	}
	newXMLNode("p",parent = top)
	div<-newXMLNode("div",attrs=c(id="show",style="text-align: center;"), "test",parent = top)

	for (nm in 2:length(labelsmut)){
					
		
		for (nr in 1:8){
			iden<-paste(sep="",nm,nr)
# 			p<-newXMLNode("p",paste(labelsmut[nm],labelsrow[nr]),parent = top)
			#tbl<-newXMLNode("table",attrs=c(id=iden,align= "center",border=1), parent = top)
			tbl<-newXMLNode("table",attrs=c(style="visibility:hidden;",id=iden,align= "center",border=1), parent = top)
			

			for (ns in 1:length(names)){
				index<-ns*8-8+nr
				temp<-data.frame(chr=(listmu[[nm]][[index]]),flag=(listfreqmu[[nm]][[index]]))
				if (nrow(temp)==0){

					temp<-data.frame(chr=0,flag=0)
				}
# 				print (temp)
				names(temp)<-c("chr",names[ns])
				
				if (ns==1){
					table<-temp
					
				}else{
					table<-merge(temp,table,by="chr",all=TRUE)
					
				}
			}
			table<-table[table$chr!=0,]
			table[is.na(table)]<-0
			tr<-newXMLNode("tr",parent = tbl)
			td<-newXMLNode("th",attrs=c(bgcolor=colall[nr]),parent=tr)
			fontcolor="black"
			if (nr>=3 & nr<=6){fontcolor="white"}
			font<-newXMLNode("font",attrs=c(color=fontcolor),labelsrow[nr],parent=td)
			col=lapply(names,function(x) newXMLNode("th", x,parent=tr))
			addChildren(tr, col)
# 			print(table)
			if (nrow(table)>0){
			for (nr in 1:nrow(table)){
				tr<-newXMLNode("tr",parent = tbl)
# 				print(table[nr,])
				text<-as.character(table[nr,1])
# 				print(text)
				color<-"white"
				if (prod(table[nr,2:ncol(table)])>0) {
					color<-"#FF6666"
				}
				td<-newXMLNode("td", attrs=c(bgcolor=color),text,parent=tr)
				text<-table[nr,2:ncol(table)]
				col=lapply(text,function(x) newXMLNode("td", x,parent=tr))
				addChildren(tr, col)

			}
			}
		}
	
	}

	
	saveXML(top,file=paste(sep="",list$path,"result.html"))


}






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

temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])

size<-3

if (is.null(list$size)==F){
	
	size<-as.numeric(list$size)
}
if (is.null(list$ini)==F){
	
	ini<-as.numeric(list$ini)
}
if (is.null(list$error)==T){
	list$error=FALSE
}
cof<-as.numeric(list$cof)/100
typemut<-0
labelsmut<-0
listmu<-vector("list")
listfreqmu<-vector("list")
nm<-2
listmu[[1]]<-0
if (is.null(list$trimmed5)==F){
	size<-size*2
	typemut<-c(typemut,4)
	labelsmut<-c(labelsmut,"5 trimming")
	listmu[[nm]]<-vector("list",length=8*length(list$group1))
	listfreqmu[[nm]]<-vector("list",length=8*length(list$group1))
	nm<-nm+1
}
if (is.null(list$trimmed3)==F){
	size<-size*2
	typemut<-c(typemut,5)
	labelsmut<-c(labelsmut,"3 trimming")
	listmu[[nm]]<-vector("list",length=8*length(list$group1))
	listfreqmu[[nm]]<-vector("list",length=8*length(list$group1))
	nm<-nm+1
}
if (is.null(list$addition3)==F){
	
	typemut<-c(typemut,6)
	labelsmut<-c(labelsmut,"3 addition ")
	listmu[[nm]]<-vector("list",length=8*length(list$group1))
	listfreqmu[[nm]]<-vector("list",length=8*length(list$group1))
	nm<-nm+1
}
if (is.null(list$mut)==F){
	
	typemut<-c(typemut,7)
	labelsmut<-c(labelsmut,"nt substitution")
	listmu[[nm]]<-vector("list",length=8*length(list$group1))
	listfreqmu[[nm]]<-vector("list",length=8*length(list$group1))
	nm<-nm+1
}
# listmu
ns<-1
listsamples<-vector("list",length=length(list$group1))
lensamples<-0
listperfect<-vector("list",length=length(list$group1))
listma<-vector("list",length=length(list$group1))
ma<-matrix(ncol=(length(typemut)-1),nrow=8)
indtypemut<-1:length(typemut)
ma[is.na(ma)]<-0
checktrimming5<-0
checktrimming3<-0
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	query<-paste(sep="","select `id`,`chr`,`freq`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`trimmed3ref`,`trimmed5ref` from `",s,"` ",options," AND `chr` NOT LIKE 'amb' AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+' AND `mut` NOT REGEXP '^[0-9]+[ACTG]+[0-9]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
# 	print(temp[temp$trimmed5ref=="na"])
	listsamples[[ns]]<-temp
	if (length(temp$id[temp$trimmed5ref!="na"])==0 & nrow(temp)>0){
		checktrimming5<-1
	}
	if (length(temp$id[temp$trimmed3ref!="na"])==0 & nrow(temp)>0){
		checktrimming3<-1
	}
	lensamples[ns]<-length(unique(listsamples[[ns]][,2]))
	listperfect[[ns]]<-temp[temp$trimmed5=="0" & temp$trimmed5!="na" & temp$trimmed3=="0" & temp$trimmed3!="na" & temp$addition3=="0" & temp$addition3!="na" & temp$mut=="0" & temp$mut!="na",2:3]
# 	print (listperfect[ns])
	if (ns==1){
		listmi<-unlist(unique(listsamples[[ns]][,2]))
	}else{
		listmi<-unlist(unique(c(listmi,listsamples[[ns]][,2])))
	}
	listma[[ns]]<-ma
	ns<-ns+1
}
# lensamples[1]
# listsamples[[1]]


# ma<-as.data.frame(ma)
# print(ma)
nts<-c("A","T","C","G")
ntspos<-c(1,2,3,4)

# print(checktrimming5)
# print(length(listmi))
if (length(listmi)>0 & length(typemut)>1){
for (chr in listmi){
#  	print (chr)
	for (ns in 1:length(list$group1)){
#  		print (ns)
		tempperfect<-listperfect[[ns]]
#  		print (tempperfect[tempperfect$chr==chr,])
		if (length(tempperfect$freq[tempperfect$chr==chr])==0){
			freqperfect<-0
		}else{
			freqperfect<-tempperfect$freq[tempperfect$chr==chr]
		}
		#print (freqperfect)
		for (nm in typemut[2:length(typemut)]){
#  			print(typemut)
			tempdata<-listsamples[[ns]]
			
			if (length(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])>0){
#  				print (tempdata)
				tempvalue<-sum(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])
				minvalue<-min(freqperfect,tempvalue,na.rm=T)
# 				print(c(freqperfect,tempvalue))
# 				temp<-tempdata
				temp<-tempdata[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na",]
# 				temp$freq<-as.numeric(temp$freq)
# 				print(list$error)
 		#		print (temp)
				temp<-filter(temp,minvalue*cof,nm,"filter",freqperfect,list$error)
# 				print (temp)
				freqperfect<-freqperfect+(tempvalue-sum(temp$freq,na.rm=T))
				if (length(temp$id)>0){
				tempma<-matrix(nrow=4,ncol=25)
				tempma[is.na(tempma)]<-0
# 				print (temp)
# 				print (nrow(temp))
#  				print (temp)
				for (i in 1:length(temp$id)){	
					if (nm<=5){
						
						mutation<-unlist(strsplit(temp[i,nm],""))
						if (mutation[1]=="q"){
							index<-13+length(mutation)-1
							tempma[1,index]<-tempma[1,index]+temp$freq[i]
							
						}else{
							index<-13-length(mutation)+1
							tempma[1,index]<-tempma[1,index]+temp$freq[i]
						}
# 						print(mutation)
						
					}else if(nm==6){
						mutation<-unlist(strsplit(temp[i,nm],""))
						tempma[1,1:(length(mutation)-1)]<-tempma[1,1:(length(mutation)-1)]+temp$freq[i]
						
					}else{
						mutation<-unlist(strsplit(gsub("[0-9]+","",temp[i,nm]),""))
						
						pos<-as.numeric(strsplit(gsub("[ATGC]+"," ",temp[i,nm])," "))
						countnt<--1
						for (posmutmir in pos){
							countnt<-countnt+2
							nt<-mutation[countnt]
							
							tempma[ntspos[nts==nt],posmutmir]<-tempma[ntspos[nts==nt],posmutmir]+temp$freq[i]
						}
 						#print(mutation)

					}
					
				}
  				#print(tempma)
				ma<-listma[[ns]]
				maxvalue<-max(tempma)
				ratio<-freqperfect/(freqperfect+maxvalue)
				#print(ratio)
				if (freqperfect+maxvalue==0){
					ratio<-1
				}
				value<-as.numeric(cut(ratio,breaks=c(-1,seq(0.2,1,0.2)),labels=1:5,right=T))
#  				print (c(freqperfect,maxvalue))
				ma[value+3,indtypemut[typemut==nm]-1]<-ma[value+3,indtypemut[typemut==nm]-1]+1

 				#print(c(value,ratio))
				index<-8*ns-8+value+3
 				templistmu<-listmu[[indtypemut[typemut==nm]]][[index]]
 				templistmu<-unlist(unique(c(templistmu,chr)))
 				listmu[[indtypemut[typemut==nm]]][[index]]<-templistmu
				listfreqmu[[indtypemut[typemut==nm]]][[index]]<-c(listfreqmu[[indtypemut[typemut==nm]]][[index]],maxvalue)
# # # # 				
				if (length(tempma[tempma>0])==1){
					ma[1,indtypemut[typemut==nm]-1]<-ma[1,indtypemut[typemut==nm]-1]+1
					index<-8*ns-8+1
					templistmu<-listmu[[indtypemut[typemut==nm]]][[index]]
					templistmu<-unlist(unique(c(templistmu,chr)))
					listmu[[indtypemut[typemut==nm]]][[index]]<-templistmu
					listfreqmu[[indtypemut[typemut==nm]]][[index]]<-c(listfreqmu[[indtypemut[typemut==nm]]][[index]],maxvalue)
				}else if (length(tempma[tempma>0])==2){
					ma[2,indtypemut[typemut==nm]-1]<-ma[2,indtypemut[typemut==nm]-1]+1
					index<-8*ns-8+2
					templistmu<-listmu[[indtypemut[typemut==nm]]][[index]]
					templistmu<-unlist(unique(c(templistmu,chr)))
					listmu[[indtypemut[typemut==nm]]][[index]]<-templistmu
					listfreqmu[[indtypemut[typemut==nm]]][[index]]<-c(listfreqmu[[indtypemut[typemut==nm]]][[index]],maxvalue)
				}else{
					ma[3,indtypemut[typemut==nm]-1]<-ma[3,indtypemut[typemut==nm]-1]+1
					index<-8*ns-8+3
					templistmu<-listmu[[indtypemut[typemut==nm]]][[index]]
					templistmu<-unlist(unique(c(templistmu,chr)))
					listmu[[indtypemut[typemut==nm]]][[index]]<-templistmu
					listfreqmu[[indtypemut[typemut==nm]]][[index]]<-c(listfreqmu[[indtypemut[typemut==nm]]][[index]],maxvalue)
				}
				listma[[ns]]<-ma
				}
				
			}
# 			print (temp)
			
		}
# 		print(ma)
# 		quit()
	}

}#for
	#print (listmu)
	writexml(listma,lensamples,list$group1,labelsmut)

}else{
	doc = newXMLDoc()
	top=newXMLNode("html")
	div<-newXMLNode("div",attrs=c(style="text-align: center;"),parent=top)
	newXMLNode("b","Check:",parent=top)
	newXMLNode("br",parent=top)
	newXMLNode("b","Remember to choose some of the checkbox options to run the analysis",parent=top)
	newXMLNode("br",parent=top)
	newXMLNode("b","If the error persits, please contact to lorena.pantano@crg.es",parent=top)
	saveXML(top,file=paste(sep="",list$path,"result.html"))

}#if error



dbDisconnect(con) 

