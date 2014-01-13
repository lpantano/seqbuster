
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
dorefseqtrm<-function(t,size){
	
	nt<-unlist(strsplit(t,""))
	s<-size/2
	nts<-length(nt)/2
	
	if (nts>=s){
		nt<-nt[(nts-s+1):(nts+s)]
# 		print(nt)
	}else{
		nt<-rep("N",size)
	}
	return(nt)
}
dorefseq<-function(seq,mut,t5,len,ini,end){
# 	print (t5)
	nt<-unlist(strsplit(seq,""))
	nt<-nt[1:len]
	if (mut!=0){
		mutation<-unlist(strsplit(gsub("[0-9]+","",mut),""))
		nt1<-mutation[1]
		nt2<-mutation[2]
		pos<-as.numeric(gsub("[ATGC]+","",mut))
		nt[pos]<-nt2
	}
	if (t5!=0){
		ntt5<-unlist(strsplit(t5,""))
		if(ntt5[1]=="q"){
			nt<-nt[(length(ntt5)-1):length(nt)]
		}else{
			nt<-c(ntt5[2:length(ntt5)],nt)
		}
	}
	if (length(nt)>=end){
		newnt<-nt[1:end]
	
	}else{
		newnt<-nt
		newnt<-append(newnt,rep("N",end-length(nt)))
		newnt<-nt[1:end]
	}
	return(newnt)
}

getmax<-function(num1,num2){
	p<-max(num1,num2)
	return (p)

}
putvalues<-function(x,y,trnt,text,n){

 tdnt<-newXMLNode("td",attrs=c(bgcolor=y),parent=trnt)
 a<-newXMLNode("a",attrs=c(href=paste(sep="",text,n)),x,parent=tdnt)

#  t<-newXMLNode("td",attrs=c(bgcolor=y), x,parent=trnt)
 
#   return(a)
}
sampleboot<-function(n,nt,totalmi,numpos,s){
	freqNTHairpinS<-freqNTHairpin[[s]]
	rn<-round(runif(totalmi)*nrow(freqNTHairpinS))
	
	rn[rn==0]<-1
	
	temp<-freqNTHairpinS[rn,]
# 	print(totalmi)
# 	print(nrow(temp))
# 	print(temp)
# 	if (nrow(temp))
	nnts<-length(temp[temp[,numpos]==nt,1])
	return(nnts)
}

getpvaluesNT<-function(numChanges,numpos,totalmi,nt,type,s){
# 	rn<-0
	p<-0
	
	if (type=="3 addition"){
# 		random<-rnorm(totalmi,mean=numChanges,sd=numChanges/10)
		changes<-pnorm(numChanges,mean=0.25*totalmi,sd=0.25*totalmi/10)
		if (changes>0.95){
				p<-1
		}
	}else{
# 		print(numpos)
# 		print(indlenreal)
		postr<-indtrlen[indlenreal==numpos]
# 		print(postr)
		
		freqNTHairpinS<-freqNTHairpin[[s]]
# 		print (c("changes",numChanges))
# 		print (c("numpos",numpos))
		if (totalmi==1){
			p<-length(freqNTHairpinS[freqNTHairpinS[,postr]==nt,1])/length(freqNTHairpinS[freqNTHairpinS[,postr]!="N",1])
		}else{
		changes<-0
		nnt<-mapply(sampleboot,1:400,nt,totalmi,postr,s)
		if (length(nnt[numChanges>=nnt])>0){
			changes<-length(nnt[numChanges>=nnt])
		}
# 		print (c("pos",numpos))
# 		print (c("nnt",nt))
# 		print (c("nc",changes/401))
# 		print (c("nmi",totalmi))
			if (changes/401>=0.975 ){
				p<-1
			}
# 		print (c("p",p))
		}
	}	
return (p)	
}

samplefunction<-function(n,numnt,nt){
# 	print (nt)
# 	print (numnt)
	pop<-1:4
	popr<-setdiff(pop,nt)
	popr[3]<-nt
	x<-sample(popr,numnt,replace =TRUE)	
	y<-length(x[x==nt])
	return (y)
}



bootfunction<-function(numnt,nt,numtotal){
	z<-0
	if (numtotal>0){
	x<-mapply(samplefunction,1:1000,numtotal,nt)
# 	print (x)
	y<-(length(x[x<=numnt]))
# 	print (y)

		if (y>0){
	# 		if (y/1001<=0.025 | y/1001>=0.975){
			if (y/1001>=0.95){
				z<-1
			}
		}else{
			z<-1
		}
	}else{
		z<-1
	}
	if (numnt==0){
		z<-0
	}
	return (z)
}

writexml<-function(listma,lensamples,names,labelsmut,labcol,numcol,pospvalue,type,ini,end){
# 	print (c(listma,lensamples,names,labelsmut,labcol,numcol,pospvalue,type,ini,end))
# 	listmaN<-listma[[1]][,numcol]

	doc = newXMLDoc()
	top=newXMLNode("body")
	div<-newXMLNode("div",attrs=c(style="text-align: center;"), parent = top)
	newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = div)

	colr<-rev(c("white","#FFCC99","#CC9966","#996633","#663300","#330000"))
	colrnames<-c("0-20","20-40","40-60","60-80","80-100")

	text<-"The variability significance is calculated using the following equation: r=Fr/(Fr+Fv), where Fr is the frequency of the reference sequence and Fv is the frequency of the variant sequence."
	newXMLNode("p",text,parent=top)
	
	tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(colr,function(x) newXMLNode("td", attrs=c(bgcolor=x,style="height:1em;"),parent=tr))
	addChildren(tr, col)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(colrnames,function(x) newXMLNode("td", x,parent=tr))
	addChildren(tr, col)
	newXMLNode("p",parent=top)
	
	div<-newXMLNode("div",attrs=c(id="show",style="text-align: center;"), "test",parent = top)
	newXMLNode("p",parent=top)
	nts<-c("A","T","C","G")
	ntspos2<-c("AA","AT","AC","AG","TA","TT","TC","TG","CA","CT","CC","CG","GA","GT","GC","GG")
	
	ntspos2ind<-1:16
	ntspos<-c(1,2,3,4)
	maxyplot<-6
	maxxplot<-2+4*length(numcol)
	coln<-c("red","blue","green","orange")
# 	rn<-runif(1)
# 	listypos<-vector("list",length(lensamples))
	xlab<-seq(3,maxxplot,4)
	ylab<-6
	
	listxpos<-seq(0.5,maxxplot-2,4)
	xposnt<-seq(1.5,maxxplot-1,1)
	xposline<-seq(1,maxxplot+1,4)
	xposntline<-seq(1,maxxplot-1,1)
	xposntlinesummary<-seq(maxxplot+3,maxxplot+7,1)
	xposntsummary<-seq(maxxplot+3.5,maxxplot+6.5,1)
# 	maxsummary<-matrix(ncol=ncol(listmaN[[1]]),nrow=length(lensamples))
# 	bitmap(paste(sep="","/srv/www/htdocs/demo2/pictures/tablemut2.bmp"),height =length(lensamples)*5, width = length(numcol)*3,units='cm',type="png16m")
	
	tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
	
	

	for (p in numcol){
		tr<-newXMLNode("tr",  parent = tbl)
		td<-newXMLNode("td", attrs=c(colspan=8,align="center",bgcolor="grey"),paste("pos",p) ,parent = tr)
		tr<-newXMLNode("tr",  parent = tbl)
		for (nt1 in nts){
		
		
# 		text<-paste(listmu[[nm]][[index]])
		text<-paste(sep="","#",p,nt1)
		td<-newXMLNode("td",attrs=c(bgcolor=coln[nts==nt1],width=2),parent = tr)

		table2<-matrix(ncol=length(names)+2,nrow=1)
		table2<-as.data.frame(table2)
		names(table2)<-c("f","chr",names)
		table2[is.na(table2)]<-"0"
		table2ratio<-table2
			for (nt2 in nts){
# 			print(nt2)
# 				for (nra in 1:5){
					for (ns in 1:length(names)){
						templab<-paste(sep="",nt1,nt2)
						index<-txs*ns-txs+indlab[templab==lab]
# 						print(templab)
# 						print(index)
						if (length(listmu[[p]][[index]])>0){
							
							temp<-data.frame(chr=listmu[[p]][[index]],freq=listfreqmu[[p]][[index]])
							tempratio<-data.frame(chr=listmu[[p]][[index]],freq=listmuratio[[p]][[index]])
							
						}else{
							temp<-data.frame(chr="0",freq=0)
							tempratio<-temp
						}
						names(temp)<-c("chr",names[ns])
						names(tempratio)<-c("chr",names[ns])
						if (ns==1){
							table<-temp
							tableratio<-tempratio
						}else{
# 							print(table)
# 							print(temp)
							table<-merge(table,temp,by="chr",all=TRUE)
							tableratio<-merge(tableratio,tempratio,by="chr",all=TRUE)
# # 							print(tableratio)
							
						}
					
					}
					table<-cbind(rep(nt2,nrow(table)),table)
					tableratio<-cbind(rep(nt2,nrow(tableratio)),tableratio)
# 					table<-cbind(rep(nra,nrow(table)),table)
# 					print (table2)
# 					print (tableratio2)
					names(table)[1]<-"f"
					names(tableratio)[1]<-"f"
					table[is.na(table)]<-0
					table2<-rbind(table2,table)
					tableratio[is.na(tableratio)]<-6
					table2ratio<-rbind(table2ratio,tableratio)
# 					print(table2)
# 				}
# 				print(table)
# 				print(table2)
# 				table2<-rbind(table2,table)
				
				#rbind table2 with r=nra
			}
# 			print(table2)
			table2<-table2[table2$chr!="0",]
			table2ratio<-table2ratio[table2ratio$chr!="0",]
			if (nrow(table2)>0){
# 			print(table2)
			td<-newXMLNode("td",parent = tr)
			a<-newXMLNode("a",attrs=c(href=text),nt1,parent=td)
			iden<-paste(sep="",p,nt1)
			tblnt<-newXMLNode("table",attrs=c(id=iden,align= "center",border=1), parent = top)
			trnt<-newXMLNode("tr",  parent = tblnt)
			text<-paste("pos",p,nt1)
			tdnt<-newXMLNode("td","NT",parent = trnt)
			col=mapply(function(x,y) newXMLNode("td", x),names)
			addChildren(trnt, col)
			tdnt<-newXMLNode("td",attrs=c(colspan=1),text,parent = trnt)
			for (i in 1:nrow(table2)){
				trnt<-newXMLNode("tr",  parent = tblnt)
				color<-as.character(colr[as.numeric(table2ratio[i,3:ncol(table2ratio)])])
				
# 				print (table2ratio[i,3:ncol(table2)])
# 				print (color)
# 				color[is.na(color)]<-"white"
				color<-c(color)
# 				print (color)
# 				print (table2[i,2:ncol(table2)])
# # 				tdnt<-newXMLNode("td",attrs=c(bgcolor=colr[as.numeric(table2[i,1])],width=2),parent = trnt)
				tdnt<-newXMLNode("td", attrs=c(bgcolor=coln[nts==table2[i,1]]),table2[i,1],parent = trnt)
				
				tdnt<-newXMLNode("td",table2[i,2],parent = trnt)
				
				
				
				mapply(putvalues,table2[i,3:ncol(table2)],color,trnt,text,names)
# 				addChildren(trnt, col)
			}
			}else{
				td<-newXMLNode("td","--",parent = tr)
			}
			
		}
	}
	
	
	jpeg(paste(sep="",list$path,"img.jpg"),height = length(lensamples)*350, width =maxxplot*80 ,units='px',quality=100)
	par(mfrow=c(length(lensamples),1),cex=1/length(names)*2,mar=c(1,1,1,1))
	for (i in 1:length(lensamples)){
		
		plot(1, type="n", axes=F,ylab="",xlim=c(0,maxxplot+10),ylim=c(0,maxyplot+3),main=names[i],xlab="")
		text(xlab,ylab,labcol)
		text(mean(xlab),maxyplot+1,"isomiR")
		text(0,4:1,c("A","U","C","G"))
		text(-0.3,3,"miRNA",srt=90,pos=2)
		text(xposnt,5,c("A","U","C","G"))
		text(xposntsummary,5,c("A","U","C","G"))
		text(mean(xposntsummary),6,"isomiR")
		text(xposntlinesummary[1],4:1,c("A","T","C","G"),pos=2)
		text(xposntlinesummary[1]-0.5,3,"miRNA",srt=90,pos=2)
		segments(xposline,5,xposline,0)
		segments(xposntline,4,xposntline,0)
		segments(xposntlinesummary,4,xposntlinesummary,0)
		totalmipos<-apply(listma[[i]][[1]],2,sum)+apply(listma[[i]][[2]],2,sum)+apply(listma[[i]][[3]],2,sum)+apply(listma[[i]][[4]],2,sum)
		totalmipos<-totalmipos[numcol]
# 		print(totalmipos)
		
		for (nt in 1:4){
			
			listma[[i]][[nt]]<-listma[[i]][[nt]][,numcol]
			ypos<-5-nt
# 			print (apply(listma[[i]][[nt]],2,sum))
# 			print(apply(listma[[i]][[nt]],2,sum))
			
			pvectref<-mapply(getpvaluesNT,apply(listma[[i]][[nt]],2,sum),labcol,totalmipos,nts[nt],type,i)
			fontv<-rep(1,length(numcol))
# 			print("ok")
			
			fontv[pvectref==1]<-2
# 			print(c("ref",fontv))
# 			print(fontv)
# 			print(pvectref)
			for (nt2 in 1:4){
				
				xpos<-listxpos+nt2
				####pvalue ref
				tfontv<-rep(1,length(numcol))
				tfontv[listma[[i]][[nt]][nt2,]>0]<-fontv[listma[[i]][[nt]][nt2,]>0]
# 				print(c("mut and red",tfontv))
				pvect<-mapply(bootfunction,listma[[i]][[nt]][nt2,],nt,apply(listma[[i]][[nt]],2,sum))
				if (length(pvect[pvect==1])>0){
					text(xpos[pvect==1],ypos,listma[[i]][[nt]][nt2,pvect==1],col="red",font=tfontv[pvect==1])
# 					text(xposntsummary[nt2],ypos,sum(listma[[i]][[nt]][nt2,]))
					text(xposntsummary[nt2],ypos,sum(listma[[i]][[nt]][nt2,pvect==1 & tfontv==2]))
					
				}
				if (length(pvect[pvect!=1])>0){
					text(xpos[pvect!=1],ypos,listma[[i]][[nt]][nt2,pvect!=1],font=tfontv[pvect!=1])
				}
			}
			
		}
		temppos<-0
# 		print(listchanges)
		
		
		
		
	}	
	dev.off()
	
	
	saveXML(top,file=paste(sep="",list$path,"result.html"))


	

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

if (is.null(list$linker)==F){
	
	options<-paste(sep="",options," `tag` LIKE '",list$linker,"' AND")
}
if (is.null(list$DB)==F){
	
	options<-paste(sep="",options," `DB` LIKE '",list$DB,"' AND")
}

if (is.null(list$locus)==F){
	
	options<-paste(sep="",options," `chr` LIKE '",list$locus,"' AND")
}

temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])



if (is.null(list$end)==F){
	
	size<-as.numeric(list$end)
}
if (is.null(list$start)==F){
	
	ini<-as.numeric(list$start)
}
if (is.null(list$error)==T){
	list$error=FALSE
}
cof<-as.numeric(list$cof)/100
# typemut<-0
#labelsmut<-0
typemut<-c(0,7)
labelsmut<-c(0,"mutation")
numcol<-c(ini:size)
labcol<-c(ini:size)
pospvalue<-c(ini:size)
iniref<-5
endref<-8
type<-"mut"

freqNTHairpin<-vector("list",length=length(list$group1))
ns<-1
listsamples<-vector("list",length=length(list$group1))
lensamples<-0
listperfect<-vector("list",length=length(list$group1))
listma<-vector("list",length=length(list$group1))
listprofile<-vector("list",length=length(list$group1))
listmitable<-vector("list",length=length(list$group1))
listchanges<-vector("list",length=16*25)
listmitableall<-listchanges

ma<-matrix(ncol=25,nrow=4)
ma[is.na(ma)]<-0
listNT<-vector("list",length=4)
tempmaO<-vector("list",length=4)
listNT[[1]]<-ma
tempmaO[[1]]<-vector("list",length=4)
listNT[[2]]<-ma
tempmaO[[2]]<-vector("list",length=4)
listNT[[3]]<-ma
tempmaO[[3]]<-vector("list",length=4)
listNT[[4]]<-ma
tempmaO[[4]]<-vector("list",length=4)
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	if (type=="mut"){
		
		query<-paste(sep="","select`chr`,`seq`,`mut`,`start`,`end`,`trimmed5` from `",s,"` ",options," AND `amb`=1 AND (`trimmed5` like '0' OR  `trimmed5` like 'na') AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+' GROUP BY `chr`;")
		
		rs <- dbSendQuery(con,query) 
		temp <- as.data.frame(fetch(rs))
# 		print((temp))
		table<-matrix(nrow=nrow(temp),ncol=size)
		table[,]<-mapply(dorefseq,temp$seq,temp$mut,temp$trimmed5,temp$end-temp$start+1,ini,size)
# 		print ("ok")
		freqNTHairpin[[ns]]<-table
		indlenreal<-1:length(freqNTHairpin[[ns]][1,])
		indtrlen<-1:length(indlenreal)
# 		print(freqNTHairpin[[ns]][1:10,])
	}
	query<-paste(sep="","select `id`,`chr`,`freq`,`trimmed5`,`trimmed3`,`addition3`,`mut` from `",s,"` ",options," AND `amb`=1  AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	listsamples[[ns]]<-temp
	lensamples[ns]<-length(unique(listsamples[[ns]][,2]))
	listperfect[[ns]]<-temp[temp$trimmed5=="0" & temp$trimmed5!="na" & temp$trimmed3=="0" & temp$trimmed3!="na" & temp$addition3=="0" & temp$addition3!="na" & temp$mut=="0" & temp$mut!="na",2:3]
# 	print (temp)
	if (ns==1){
		listmi<-unique(listsamples[[ns]][,2])
	}else{
		listmi<-unique(c(listmi,listsamples[[ns]][,2]))
	}
	listmitable[[ns]]<-listchanges
	listma[[ns]]<-listNT
	ns<-ns+1
}
# lensamples[1]
# listsamples[[1]]
nts<-c("A","T","C","G")
ntspos2<-c("AA","AT","AC","AG","TA","TT","TC","TG","CA","CT","CC","CG","GA","GT","GC","GG")
ntspos2ind<-1:16
ntspos<-c(1,2,3,4)
listmu<-vector("list",length=25)
listfreqmu<-vector("list",length=25)
listmuratio<-vector("list",length=25)

templistfreq2<-vector("list",length=25)
lab<-0
for (nt in ntspos2){

# for (i in 1:5){

lab<-c(lab,paste(sep="",nt))
# }

}
lab
lab<-lab[2:length(lab)]
indlab<-1:length(lab)
txs<-length(lab)
for (p in numcol){
# 	print (p)
	listmu[[p]]<-vector("list",length=txs*length(list$group1))
	listmuratio[[p]]<-vector("list",length=txs*length(list$group1))
	listfreqmu[[p]]<-vector("list",length=txs*length(list$group1))
	templistfreq2[[p]]<-vector("list",length=length(ntspos2ind))
}

for (chr in listmi){
# 	print (chr)
	for (ns in 1:length(list$group1)){
# 		print (ns)
		
		tempperfect<-listperfect[[ns]]
# 		print (tempperfect)
		if (length(tempperfect$freq[tempperfect$chr==chr])==0){
			freqperfect<-0
		}else{
			freqperfect<-tempperfect$freq[tempperfect$chr==chr]
		}
		for (nm in typemut[2:length(typemut)]){
# 			print(nm)
			tempdata<-listsamples[[ns]]
			if (length(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])>0){
				
				tempvalue<-sum(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])
				minvalue<-min(freqperfect,tempvalue,na.rm=T)
# 				print(c(freqperfect,tempvalue))
				temp<-tempdata[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na",]
				temp<-filter(temp,minvalue*cof,nm,"filter",freqperfect,list$error)
# 				freqperfect<-freqperfect+(tempvalue-sum(temp$freq,na.rm=T))
				if (length(temp$id)>0){
				
				tempma<-tempmaO
# 				templistfreq<-tempmaO
				templistfreq<-templistfreq2
				for (i in 1:length(temp$id)){	
				
						mutation<-unlist(strsplit(gsub("[0-9]+","",temp[i,nm]),""))
						pos<-as.numeric(strsplit(gsub("[ATGC]+"," ",temp[i,nm])," "))
# 						print(temp[i,nm])
						if (pos>=ini & pos<=size){
						countnt<--1
						for (posmutmir in pos){
							countnt<-countnt+2
							nt1<-mutation[countnt]
							nt2<-mutation[countnt+1]
							tempma[[ntspos[nts==nt2]]][[ntspos[nts==nt1]]]<-unique(append(tempma[[ntspos[nts==nt2]]][[ntspos[nts==nt1]]],posmutmir))
							templab<-paste(sep="",nt2,nt1)
# 						print(tempma)
							templistfreq[[posmutmir]][[ntspos2ind[ntspos2==templab]]]<-sum(templistfreq[[posmutmir]][[ntspos2ind[ntspos2==templab]]],temp$freq[i])
# 						if (chr=="hsa-miR-23b"){
# 							print(temp[i,nm])
# 							print(temp$freq[i])}
						}
						}
						# 						tempma[ntspos[nts==nt2]+4,pos]<-tempma[ntspos[nts==nt2]+4,pos]+temp$freq[i]
	
				}
				ma<-listma[[ns]]
# 				print (chr)
				for (nt1 in 1:4){
					for (nt2 in 1:4){
					if (length(tempma[[nt1]][[nt2]])>0){
# 					print(tempma[[nt1]][[nt2]])
						ma[[nt1]][nt2,tempma[[nt1]][[nt2]]]<-ma[[nt1]][nt2,tempma[[nt1]][[nt2]]]+1
						for (p in (tempma[[nt1]][[nt2]])){
# 							p<-as.numeric(pos)
							templab<-paste(sep="",nts[nt1],nts[nt2])
# 							print(p)
							ratio<-freqperfect/(freqperfect+templistfreq[[p]][[ntspos2ind[ntspos2==templab]]])
							if (is.na(ratio)==TRUE){
								ratio<-1
							}
# 							print(templistfreq[[pos]][[ntspos2ind[ntspos2==templab]]])
							value<-as.numeric(cut(ratio,breaks=c(-1,seq(0.2,1,0.2)),labels=1:5,right=T))
# 							print(c(value,ratio,freqperfect,templistfreq[[p]][[ntspos2ind[ntspos2==templab]]]))
							templab2<-paste(sep="",nts[nt1],nts[nt2])
							index<-txs*ns-txs+indlab[templab2==lab]
							listmu[[p]][[index]]<-c(listmu[[p]][[index]],chr)
# 							print(index)
							listmuratio[[p]][[index]]<-c(listmuratio[[p]][[index]],value)
							listfreqmu[[p]][[index]]<-c(listfreqmu[[p]][[index]],templistfreq[[p]][[ntspos2ind[ntspos2==templab]]])
						}
					}
					}
				
				}
				listma[[ns]]<-ma
				}
# 				quit()
			}
# 			print (temp)
			
		}
# 		print(ma)
# 		quit()
	}

}
# listfreqmu
# listmu
# listmuratio
# q()
writexml(listma,lensamples,list$group1,labelsmut,labcol,numcol,pospvalue,type,iniref,endref)


dbDisconnect(con) 
