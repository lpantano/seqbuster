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

getmax<-function(num1,num2){
	p<-max(num1,num2)
	return (p)

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
# 		changes<-0
		nnt<-mapply(sampleboot,1:400,nt,totalmi,postr,s)
		if (length(nnt[numChanges>=nnt])>0){
			changes<-length(nnt[numChanges>=nnt])
		}
# 		print (c("nnt",nnt))
# 		print (c("nc",numChanges))
# 		print (c("nmi",totalmi))
			if (changes/401>=0.95){
				p<-1
			}
			
		}
	}	
return (p)	
}

# mapply(getpvaluesNT,listma[[1]][1,],8:13,apply(listma[[1]][5:9,],2,sum),'A')


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
writexml<-function(listma,listmifinal,names,labelsmut,labcol,numcol,pospvalue,type,iniref,endref,sortt){
# 	print (numcol)
	maxyplot<-length(listmifinal)*2.2+3
	maxxplot<-length(numcol)*1.2+1
	labnt<-c("A","U","C","G")
	coln<-c("red","blue","green","yellow","red","blue","green","yellow")
	colr<-rev(c("white","#FFCC99","#CC9966","#996633","#663300","#330000"))
	colall<-c(coln,colr)
# 	if (type=="mut"){colr<-coln}
	widthpic<-length(numcol)*5
# 	bitmap(paste(sep="","/srv/www/htdocs/demo2/pictures/",type,".bmp"),height = 30, width = widthpic,units='cm',type="png16m")
	factor<-4
	if (length(numcol)>6){
		factor<-7
	
	}
	doc = newXMLDoc()
	top=newXMLNode("html")
	newXMLNode("h3",attrs=c(align="center"),paste("Studying the variability of ",labelsmut[2]),parent=top)
			
	tbl<-newXMLNode("table",attrs=c(width="200",align="center",border=1), parent = top)
	tr<-newXMLNode("tr",  parent = tbl)
	td<-newXMLNode("td", attrs=c(colspan=2,bgcolor="orange",align="center"),"Cell scheme", parent = tr)
	tr<-newXMLNode("tr",  parent = tbl)
	td<-newXMLNode("td", "Nucleotide code", parent = tr)
	td<-newXMLNode("td",  "Variability significance code",parent = tr)
	newXMLNode("p",parent=top)
	colrnames<-c(">80","60-80","60-40","40-20","<20")

	br<-newXMLNode("p","Percentage of the isomiRs with respect to the corresponding reference miRNAs. ",attrs=c(align="center"),parent=top)
	
	text<-"The variability significance is calculated using the following equation: r=Fv/(Fr+Fv)*100, where Fr is the frequency of the reference sequence and Fv is the frequency of the variant sequence."
	newXMLNode("a","?",attrs=c(title=text,href=""),parent=br)
	tbl<-newXMLNode("table",attrs=c(align="center",border=1), parent = top)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(colall[9:13],function(x) newXMLNode("td", attrs=c(bgcolor=x,style="height:1em;"),parent=tr))
	addChildren(tr, col)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(colrnames,function(x) newXMLNode("td", x,parent=tr))
	addChildren(tr, col)
	
	div<-newXMLNode("div",attrs=c(style="text-align: center;"),parent=top)
	newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = div)
# 	top=newXMLNode("div",attrs=c(style="text-align: left;"))
	
	jpeg(paste(sep="",list$path,"img.jpg"),height = maxyplot*15+3*15, width =widthpic*50 ,units='px',quality = 100)
	par(oma=c(2,1,3,1),mar=c(3,1,7,1),plt=c(0.03,0.99,0.03,0.99))
	plot(1, type="n", axes=F, xlab="", ylab="",xlim=c(-1,maxxplot+1),ylim=c(0,maxyplot+3),main="",xaxs="i",yaxs="i")
	legend("top",labnt,fill=coln, cex=1/length(numcol)*factor,horiz=TRUE)
	xvar<-(2)/5
	poslegend<-seq(maxxplot/2-1,maxxplot/2+1,xvar)
	yvar<-((maxyplot+3)-(maxyplot+1.5))/5
# 	print(poslegend)
	ylegend<-seq((maxyplot),(maxyplot+1.5),yvar)
# 	print(ylegend)
# 	if (iniref!=5){
		polygon(c(poslegend[1],poslegend[2],poslegend[2]),c(ylegend[1],ylegend[1],ylegend[2]),col=colr[5],border='NA')
		polygon(c(poslegend[2],poslegend[2],poslegend[3],poslegend[3]),c(ylegend[1],ylegend[2],ylegend[3],ylegend[1]),col=colr[4],border='NA')
		polygon(c(poslegend[3],poslegend[3],poslegend[4],poslegend[4]),c(ylegend[1],ylegend[3],ylegend[4],ylegend[1]),col=colr[3],border='NA')
		polygon(c(poslegend[4],poslegend[4],poslegend[5],poslegend[5]),c(ylegend[1],ylegend[4],ylegend[5],ylegend[1]),col=colr[2],border='NA')
		polygon(c(poslegend[5],poslegend[5],poslegend[6],poslegend[6]),c(ylegend[1],ylegend[5],ylegend[6],ylegend[1]),col=colr[1],border='NA')
		text(mean(poslegend),-3,"editting significance",pos=1,cex=1/length(numcol)*factor)
# 		#group of mutation
# 	}
 	ypos<-seq(0,maxyplot,2.2)
 	xpos<-seq(1,maxxplot-1,1.2)

	script<-newXMLNode("script", attrs=c(language="JavaScript",src="showtableid.js"),"",parent=top)
	newXMLNode("p",attrs=c(align="center"),"Cell colors indicate the abundance of isomiRs (%) with respect to the corresponding reference miRNA",parent=top)
	newXMLNode("p",parent=top)
	div<-newXMLNode("div",attrs=c(id="show",style="text-align: center;"), "Click on name to show further information",parent = top)
	newXMLNode("p",parent=top)	
	tblbig<-newXMLNode("table",attrs=c(align= "center",border=0), parent = top)
	
	trbig<-newXMLNode("tr",  parent = tblbig)
	
	
	indrows<-1
	ipos<-0
	text(xpos+0.5,ypos[length(ypos)-1],labcol,pos=3)
	for (i in sortt$ix){
# 	for (i in sortt$ix[1:4]){
		if (indrows==11){
			indrows<-1
			trbig<-newXMLNode("tr",  parent = tblbig)
			
		}
		indrows<-1+indrows
		#text<-paste(sep="","#",i)
        text<-paste(sep="","javascript:loadtable('",i,"')")
		tdbig<-newXMLNode("td", attrs=c(colspan=length(labcol)), parent = trbig)
		a<-newXMLNode("a",attrs=c(href=text),listmifinal[i],parent=tdbig)
		ipos<-ipos+1
# 		print(listmifinal[i])
# 		r<-2
# 		if(list$nt=="chang"){
# 			r<-6
# 		}
		text(1,(ypos[ipos]+1),listmifinal[i],pos=2,cex=1/length(numcol)*factor)
		tempy1<-rep(ypos[ipos],length(xpos))
# 		tempy2<-rep(ypos[i+1],length(xpos))
		
		#newXMLNode("p",parent = top)
		iden<-i
        tbl<-newXMLNode("table",attrs=c(style="visibility:hidden;",id=iden,width="20%",align= "center",border=1), parent = top)
		#tbl<-newXMLNode("table",attrs=c(id=iden,width="20%",align= "center",border=1), parent = top)
		tr<-newXMLNode("tr",  parent = tbl)
# 		text<-paste(listmi[i],listperfectfreq[i])
		td<-newXMLNode("td", attrs=c(bgcolor="orange"),listmifinal[i], parent = tr)
		td<-newXMLNode("td", attrs=c(bgcolor="white"),as.character(listmifinalfreq[i]), parent = tr)
		#norm listma[[i]]
		ma<-listma[[i]][1:4,numcol]
# 		maref<-listma[[i]][(5):(8),numcol]
		maA<-ma
 		ma<-scale(ma,center=FALSE,scale=colSums(ma))*2
		#ma<-ma/sum(ma)*2
		ma[is.na(ma)]<-0
		mar<-listma[[i]][9:12,numcol]
# 		print(listma[[i]][,numcol])
		mar[mar==0]<-6
		for (c in 1:ncol(ma)){
			td<-newXMLNode("td", attrs=c(bgcolor="orange"),labcol[c], parent = tr)
		}
# 		print(mar)
		tempy2<-tempy1+ma[1,]
# 		
		rect(xpos,tempy1,xpos+1,tempy1+2)
		if (max(tempy2)>max(tempy1)){
# 			print("ok")
			rect(xpos[tempy2>tempy1],tempy1[tempy2>tempy1],xpos[tempy2>tempy1]+0.5,tempy2[tempy2>tempy1],col=coln[1],border = 'NA')
# 			print(ma[1,])
			tempcr<-mar[1,tempy2>tempy1]
# 			print(mar[1,])
			rect(xpos[tempy2>tempy1]+0.5,tempy1[tempy2>tempy1],xpos[tempy2>tempy1]+1,tempy2[tempy2>tempy1],col=colr[tempcr],border = 'NA')
			
			
		}
		
		tr<-newXMLNode("tr",  parent = tbl)
		td<-newXMLNode("td",attrs=c(bgcolor="yellow",colspan=2),labnt[1], parent = tr)
		for (c in 1:ncol(ma))
		{
			tempcr<-mar[1,c]
# 			print (tempcr)
			colorfont<-"black"
			if (tempcr<=3){
				colorfont<-"white"
			}
			td<-newXMLNode("td",attrs=c(bgcolor=colr[tempcr]),parent = tr)
			font<-newXMLNode("font",attrs=c(color=colorfont),as.character(round(maA[1,c],digits=2)),parent=td)
			
			
		}
		counter<-2
		for (j in 2:4){
			tempy1<-tempy2
			tempy2<-tempy1+(ma[j,])
			tr<-newXMLNode("tr", parent = tbl)
			td<-newXMLNode("td",attrs=c(bgcolor="yellow",colspan=2),labnt[j], parent = tr)
			for (c in 1:ncol(ma))
			{
				tempcr<-mar[counter,c]
# 				print(tempcr)
				colorfont<-"black"
				
				if (tempcr<=3){
					colorfont<-"white"
				}
				
				td<-newXMLNode("td",attrs=c(bgcolor=colr[tempcr]),parent = tr)
				font<-newXMLNode("font",attrs=c(color=colorfont),as.character(round(maA[counter,c],digits=2)),parent=td)
				
			}
			
			if (max(tempy2)>max(tempy1)){
				rect(xpos[tempy2>tempy1],tempy1[tempy2>tempy1],xpos[tempy2>tempy1]+0.5,tempy2[tempy2>tempy1],col=coln[j],border ='NA')
			}
			
			if (max(tempy2)>max(tempy1)){
				tempcr<-mar[counter,tempy2>tempy1]
				rect(xpos[tempy2>tempy1]+0.5,tempy1[tempy2>tempy1],xpos[tempy2>tempy1]+1,tempy2[tempy2>tempy1],col=colr[tempcr],border = 'NA')
			
			}
# 			print(tempcr)
			counter<-counter+1
			
		}
		
		#print (names(table[c]))
	}
		
# 		}
			
	
	
	dev.off()
	
	
	saveXML(top,file=paste(sep="",list$path,"result.html"))


}




if (is.null(list$freq1)==F){
	temp<-unlist(strsplit(list$freq1," "))
	options<-paste(sep=" ","where `freq` >",list$freq1," AND")
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

temp<-unlist(strsplit(options," "))
options<-paste(collapse=" ",temp[1:(length(temp)-1)])



if (is.null(list$size)==F){
	
	size<-as.numeric(list$size)
}
if (is.null(list$start)==F){
	
	ini<-as.numeric(list$start)
}
if (is.null(list$error)==T){
	list$error=FALSE
}
cof<-as.numeric(list$cof)/100
sortt<-list$sort

# coFreq

typemut<-0
labelsmut<-0
iniref<-9
endref<-12
type<-0
if (list$isoseq=="5 trimming"){
	size<-size*2
	typemut<-c(typemut,4)
	labelsmut<-c(labelsmut,"5 trimming") 
	numcol<-c((13-size/2):12,14:(13+size/2))
	labcol<-c(-(size/2):-1,1:(size/2))
	pospvalue<-c((11-size/2):10,11:(11+size/2-1))
	type<-"trimmed5"
}
if (list$isoseq=="3 trimming"){
	size<-size*2
	typemut<-c(typemut,5)
	labelsmut<-c(labelsmut,"3 trimming")
	numcol<-c((13-size/2):12,14:(13+size/2))
	labcol<-c(-(size/2):-1,1:(size/2))
	pospvalue<-c((11-size/2):10,11:(11+size/2-1))
	type<-"trimmed3"
}	
if (list$isoseq=="3 addition"){
	size<-as.numeric(list$sizead)
	typemut<-c(typemut,6)
	labelsmut<-c(labelsmut,"3 addition")
	numcol<-c(1:size)
	labcol<-c(1:size)
	pospvalue<-0
	type<-"addition3"
}
if (list$isoseq=="nt substitution"){
	size<-as.numeric(list$end)
	typemut<-c(typemut,7)
	labelsmut<-c(labelsmut,"nt-substitution")
	numcol<-c(ini:size)
	labcol<-c(ini:size)
	pospvalue<-c(ini:size)
# 	iniref<-5
# 	endref<-8
	type<-"mut"
}

ns<-1
listsamples<-vector("list",length=length(list$group1))
lensamples<-0
listprofile<-vector("list",length=length(list$group1))
listperfect<-vector("list",length=length(list$group1))

ma<-matrix(ncol=25,nrow=12)
ma[is.na(ma)]<-0
checktrimming5<-0
checktrimming3<-0
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	query<-paste(sep="","select `id`,`chr`,`freq`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`trimmed5ref`,`trimmed3ref` from `",s,"` ",options," AND `amb`=1  AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+' AND `mut` NOT REGEXP '^[0-9]+[ACTG]+[0-9]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	if (length(temp$chr[temp$trimmed5ref!="na"])==0 & nrow(temp)>0 & list$isoseq=="trimmed5") {
		checktrimming5<-1
	}
	if (length(temp$chr[temp$trimmed3ref!="na"])==0 & nrow(temp)>0 & list$isoseq=="trimmed3"){
		checktrimming3<-1
	}
	listsamples[[ns]]<-temp
# 	print(nrow(temp))
	lensamples[ns]<-length(unique(listsamples[[ns]][,2]))
	listperfect[[ns]]<-temp[temp$trimmed5=="0" & temp$trimmed5!="na" & temp$trimmed3=="0" & temp$trimmed3!="na" & temp$addition3=="0" & temp$addition3!="na" & temp$mut=="0" & temp$mut!="na",2:3]
# 	print (listperfect[ns])
	if (ns==1){
		listmi<-unique(listsamples[[ns]][,2])
	}else{
		listmi<-unique(c(listmi,listsamples[[ns]][,2]))
	}
	listprofile[[ns]]<-0
	###########################WARNINGGGGGGGGGGG
# 	listprofile[[ns]]<-freqNTHairpin[freqNTHairpin[1,]==unique(listsamples[[ns]][,2]),2:ncol(freqNTHairpi)]
	
	ns<-ns+1
}
# lensamples[1]
# listsamples[[1]]


# ma<-as.data.frame(ma)
# print(length(listmi))
nts<-c("A","T","C","G")
ntspos<-c(1,2,3,4)
tempmaind<-1:25
listma<-vector("list")
listratio<-vector()
indmi<-1:length(listmi)
listmifinal<-vector()
listmifinalfreq<-vector()
listperfectfreq<-vector(length=length(indmi))
indlistma<-1
# length(listma)


for (chr in listmi){
# 	print (chr)
# 	listma[[indmi[listmi==chr]]]<-ma
	for (ns in 1:length(list$group1)){
# 		print (ns)
		
		tempperfect<-listperfect[[ns]]
# 		print (tempperfect)
		if (length(tempperfect$freq[tempperfect$chr==chr])==0){
			freqperfect<-0
		}else{
			freqperfect<-tempperfect$freq[tempperfect$chr==chr]
		}
		listperfectfreq[indmi[listmi==chr]]<-freqperfect
		for (nm in typemut[2:length(typemut)]){
# 			print(nm)
			tempdata<-listsamples[[ns]]
			tempma<-matrix(nrow=12,ncol=25)
			tempma[is.na(tempma)]<-0
			tempma[9:12,]<-6
# 			print(tempma)
			if (length(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])>0){
				
				tempvalue<-sum(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])
				minvalue<-min(freqperfect,tempvalue,na.rm=T)
# 				print(c(freqperfect,tempvalue))
				temp<-tempdata[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na",]
				temp<-filter(temp,minvalue*cof,nm,type="filter",freq=freqperfect,error=list$error)
# 				print(c(tempvalue,sum(temp$freq,na.rm=T)))
# 				freqperfect<-freqperfect+(tempvalue-sum(temp$freq,na.rm=T))
				if (length(temp$id)>0){
				
				
				for (i in 1:length(temp$id)){	
					if (nm==4){
						mutation<-unlist(strsplit(temp[i,nm],""))
						nt<-mutation[length(mutation)]
						if (mutation[1]=="q"){
							index<-13-length(mutation)+1
							tempma[ntspos[nts==nt],index]<-tempma[ntspos[nts==nt],index]+temp$freq[i]
						}else{
							index<-13+length(mutation)-1
							tempma[ntspos[nts==nt],index]<-tempma[ntspos[nts==nt],index]+temp$freq[i]
						}
					}else if (nm==5){
						mutation<-unlist(strsplit(temp[i,nm],""))
						nt<-mutation[length(mutation)]
# 						print(mutation)
						if (mutation[1]=="q"){
							index<-13+length(mutation)-1
							tempma[ntspos[nts==nt],index]<-tempma[ntspos[nts==nt],index]+temp$freq[i]
						}else{
							index<-13-length(mutation)+1
							tempma[ntspos[nts==nt],index]<-tempma[ntspos[nts==nt],index]+temp$freq[i]
						}
					}else if(nm==6){
					#change
						mutation<-unlist(strsplit(temp[i,nm],""))
						postemp<-1
						for (nt in mutation[2:length(mutation)]){
							tempma[ntspos[nts==nt],postemp]<-tempma[ntspos[nts==nt],postemp]+temp$freq[i]
							postemp<-postemp+1
						
						}
						
					}else{
# 						if (chr=="hsa-miR-1827"){
# 								print(temp[i,nm])
# 						}
						mutation<-unlist(strsplit(gsub("[0-9]+","",temp[i,nm]),""))
						
						pos<-as.numeric(strsplit(gsub("[ATGC]+"," ",temp[i,nm])," "))
# 						print(mutation)
						if (pos>=ini & pos<=size){
						r<-0
						r2<-4
						if(list$nt=="chang"){
							r<-4
							r2<-0
						}
						countnt<--1
						for (posmutmir in pos){
							countnt<-countnt+2
							nt1<-mutation[countnt]
							nt2<-mutation[countnt+1]
							tempma[ntspos[nts==nt2]+r,posmutmir]<-tempma[ntspos[nts==nt2],pos]+temp$freq[i]
							tempma[ntspos[nts==nt1]+r2,posmutmir]<-tempma[ntspos[nts==nt1]+4,pos]+temp$freq[i]
						}
						}
					}
					
				}
# 				ma<-tempma
# 				print(chr)
				#listma[[listmi==chr]]<-tempma
# 				maxvalue<-apply(tempma,2,max)
# 				if(chr=="hsa-miR-101"){
				for (i in 1:4){
				
					ratio<-freqperfect/(freqperfect+tempma[i,])
# 					if (chr=="hsa-miR-1827"){
# 						print(ratio)
# 						print(tempma[i,1:8])
# 						print(value)

# 					}
					if (is.na(ratio)==TRUE){
						ratio[is.na(ratio)]<-1
					}
					
# 					print (ratio[numcol])
					value<-as.numeric(cut(ratio,breaks=c(-1,seq(0.2,1,0.2)),labels=1:5,right=T))
					
					value[is.na(value)]<-6
					value[tempma[i,]==0]<-6
					
					tempma[i+8,]<-value
# 					print (c(value[numcol],ratio[numcol]))
					
				}
				
# 				print (c(chr,min(tempma[9:12,numcol])))
				
# 				}
				
				if (sum(tempma[1:4,numcol])>0){
# 					print(tempma[,numcol])
					listma[[indlistma]]<-tempma
# 					if (chr=="hsa-miR-1827"){
# 						print (tempma)
# 					}
					listmifinal[indlistma]<-chr
					listmifinalfreq[indlistma]<-freqperfect
					indlistma<-indlistma+1
					listratio<-append(listratio,min(tempma[9:12,numcol]))
				}
				
				}
# 				if (min(tempma[9:12,numcol])==0){
# 				print(tempma[1:4,numcol])
# 				}
# 				quit()
			}
			
			
# 			print (temp)
			
		}
# 		print(ma)
# 		quit()
	}

}
# listratio
if (sortt=="ratio"){
	listratiosort<-sort(listratio,index.return=TRUE)
}else{
	listratiosort<-sort(listmifinalfreq,index.return=TRUE)
}
# length(listma)
# listratiosort
# listma
# print(listmifinal)
# listmifinal[1:4]
writexml(listma,listmifinal,list$group1,labelsmut,labcol,numcol,pospvalue,type,iniref,endref,listratiosort)

dbDisconnect(con) 


