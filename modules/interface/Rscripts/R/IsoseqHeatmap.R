
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
	if (t5!=0){
		ntt5<-unlist(strsplit(t5,""))
		if(ntt5[1]=="q"){
			nt<-nt[(length(ntt5)-1):length(nt)]
		}else{
			nt<-c(ntt5[2:length(ntt5)],nt)
		}
	}
	if (mut!=0){
		mutation<-unlist(strsplit(gsub("[0-9]+","",mut),""))
		nt1<-mutation[1]
		nt2<-mutation[2]
		pos<-as.numeric(gsub("[ATGC]+","",mut))
		nt[pos]<-nt2
	}
	
	#fix position
	if (length(nt)>=end){
		newnt<-nt[1:end]
	
	}else{
		newnt<-nt
		newnt<-append(newnt,rep("N",end-length(nt)))
		newnt<-newnt[1:end]
	}
# 	print (newnt)
	return(newnt)
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
# 	print(nnts)
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
# 		print (c("changes",numChanges))
# 		print (c("numpos",numpos))
# 		print (c("nummi",totalmi))
		postr<-indtrlen[indlenreal==numpos]
# 		print(postr)
		
		freqNTHairpinS<-freqNTHairpin[[s]]
		
		if (totalmi==1){
# 			print("1")
			p<-length(freqNTHairpinS[freqNTHairpinS[,postr]==nt,1])/length(freqNTHairpinS[freqNTHairpinS[,postr]!="N",1])
		}else{
			changes<-0
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
dotds<-function(x,y,tr) {
	c<-rev(c("#996633","#663300","#330000"))
	colorfont<-"black"
	if (length(c[c==y])>0){
		colorfont<-"white"
	}
	td<-newXMLNode("td",attrs=c(bgcolor=y), parent=tr)
	font<-newXMLNode("font",attrs=c(color=colorfont),x,parent=td)
	return(td)
}
writexml<-function(listma,lensamples,names,labelsmut,labcol,numcol,pospvalue,type,iniref,endref){
	listmaN<-listma[[1]][,numcol]
	if (type!="mut"){numrow<-c(1:4)}else{numrow<-c(1:8)}
	
# 	rn<-runif(1)
	listmaN<-vector("list",length(lensamples))
	listypos<-vector("list",length(lensamples))
	listxpos<-vector("list",length(lensamples))
	summary<-matrix(ncol=length(numcol),nrow=length(lensamples))
	summaryN<-matrix(ncol=length(numcol),nrow=length(lensamples))
	summaryNminus<-matrix(ncol=length(numcol),nrow=length(lensamples))
# 	maxsummary<-matrix(ncol=ncol(listmaN[[1]]),nrow=length(lensamples))
	for (i in 1:length(lensamples)){
		listma[[i]]<-listma[[i]][,numcol]
# 		print(listma[[i]])
# 		print(normalize[[i]][numcol])
# 		normalize[[i]][numcol]<-rep(sum(listma[[i]][1:4,]),6)
		listmaN[[i]]<-scale(listma[[i]],scale=normalize[[i]][numcol],center=F)
		listmaN[[i]]<-listmaN[[i]]*100
# 		print(listmaN[[i]])
# 		listmaN[[i]]<-listma[[i]]/lensamples[i]*100
		listxpos[[i]]<-seq(i,i+(length(lensamples)+2)*(ncol(listmaN[[1]]))-1,length(lensamples)+2)
		listypos[[i]]<-listmaN[[i]]
		listypos[[i]][1,]<-listmaN[[i]][1,]
		listypos[[i]][2,]<-apply(listmaN[[i]][1:2,],2,sum)
		listypos[[i]][3,]<-apply(listmaN[[i]][1:3,],2,sum)
		listypos[[i]][4,]<-apply(listmaN[[i]][1:4,],2,sum)
		listypos[[i]][5,]<-listmaN[[i]][5,]
		listypos[[i]][6,]<-apply(listmaN[[i]][5:6,],2,sum)
		listypos[[i]][7,]<-apply(listmaN[[i]][5:7,],2,sum)
		listypos[[i]][8,]<-apply(listmaN[[i]][5:8,],2,sum)
		listypos[[i]][9,]<-listmaN[[i]][9,]
		listypos[[i]][10,]<-apply(listmaN[[i]][9:10,],2,sum)
		listypos[[i]][11,]<-apply(listmaN[[i]][9:11,],2,sum)
		listypos[[i]][12,]<-apply(listmaN[[i]][9:12,],2,sum)
		listypos[[i]][13,]<-apply(listmaN[[i]][9:13,],2,sum)
		
		summary[i,]<-apply(listma[[i]][1:4,],2,sum)
		summaryN[i,]<-apply(listmaN[[i]][1:4,],2,sum)
		summaryNminus[i,]<-apply(listmaN[[i]][9:13,],2,sum)
		if (type=="mut"){
			summaryNminus[i,]<-apply(listmaN[[i]][5:8,],2,sum)
		}
		
# 		summaryN[i,]<-apply(listmaN[[i]][1:4,],2,sum)
		
	}	
	factor<-4
	if (length(numcol)>6){
		factor<-9
	
	}
	maxyplot<-max(summaryN)+max(summaryN)*0.20
	maxyplotminus<-max(summaryNminus)+max(summaryNminus)*0.20
	if (max(summaryN)<10){
		maxyplot<-12
		maxyplotminus<-12
		factor<-5
	}
# 	print(summaryN)
	
	maxxplot<-(length(lensamples)+2)*ncol(listma[[1]])+1
# 	print(maxxplot)
# 	print(maxyplotminus)
# 	print(maxyplot)
	coln<-c("red","blue","green","yellow")
	colr<-rev(c("#FFCC99","#CC9966","#996633","#663300","#330000"))
	colall<-c(coln,coln,colr)
	if (type=="mut"){colr<-coln}
	widthpic<-length(numcol)/2*12
# 	bitmap(paste(sep="","/srv/www/htdocs/demo2/pictures/",type,".bmp"),height = 30, width = widthpic,units='cm',type="png16m")
	jpeg(paste(sep="",list$path,"img.jpg"),height = 500, width =widthpic*15 ,units='px',quality = 100)
	par(mar=rep(0, 4))
	plot(1, type="n", axes=F, xlab="", ylab="",xlim=c(-1,maxxplot+1),ylim=c(-maxyplotminus-30,maxyplot+0.30*maxyplot),main="")
	legend("top",c("A","U","C","G"),fill=coln, cex=1/length(numcol)*factor,horiz=TRUE)
	#xaxis
	segments(0,0,maxxplot,0)
	#yaxis
	segments(0,maxyplot,0,0)
	if (maxyplot>10) {
		yaxes<-seq(10,maxyplot,10)
		yaxesminus<-seq(-10,-maxyplotminus,-10)
	}else{
		yaxes<-10
		yaxesminus<--10
	}
	segments(0,yaxes,-0.2,yaxes)
	text(maxxplot,maxyplot/2,"% of miRNAs",srt=90,cex=1/length(numcol)*factor)
	text(-0.2,yaxes,labels=yaxes,pos=2,cex=1/length(numcol)*factor)
	#yaxisminus
	segments(maxxplot,-maxyplotminus,maxxplot,0)
	
	segments(maxxplot,yaxesminus,maxxplot+0.2,yaxesminus)
	text(0,-maxyplot/2,"% of miRNAs",srt=90,cex=1/length(numcol)*factor)
	text(maxxplot+0.2,yaxesminus,labels=-yaxesminus,pos=4,cex=1/length(numcol)*factor)
	#legend
	if (length(numcol)<4){
		xvar<-((maxxplot/2+0.2*maxxplot)-(maxxplot/2-0.2*maxxplot))/5
		poslegend<-seq(maxxplot/2-0.25*maxxplot,maxxplot/2+0.25*maxxplot,xvar)
	}else{
		xvar<-(3.5)/5
		poslegend<-seq(maxxplot/2-1,maxxplot/2+4.5,xvar)
	}
	yvar<-((maxyplotminus+15)-(maxyplotminus+20))/5
# 	print(yvar)
	ylegend<-seq(-maxyplotminus-15,-maxyplotminus-20,yvar)
# 	print(ylegend)
	if (iniref!=5){
		polygon(c(poslegend[1],poslegend[2],poslegend[2]),c(ylegend[1],ylegend[1],ylegend[2]),col=colr[5],border='NA')
		polygon(c(poslegend[2],poslegend[2],poslegend[3],poslegend[3]),c(ylegend[1],ylegend[2],ylegend[3],ylegend[1]),col=colr[4],border='NA')
		polygon(c(poslegend[3],poslegend[3],poslegend[4],poslegend[4]),c(ylegend[1],ylegend[3],ylegend[4],ylegend[1]),col=colr[3],border='NA')
		polygon(c(poslegend[4],poslegend[4],poslegend[5],poslegend[5]),c(ylegend[1],ylegend[4],ylegend[5],ylegend[1]),col=colr[2],border='NA')
		polygon(c(poslegend[5],poslegend[5],poslegend[6],poslegend[6]),c(ylegend[1],ylegend[5],ylegend[6],ylegend[1]),col=colr[1],border='NA')
		text(mean(poslegend),-20-maxyplotminus,"variability significance",pos=1,cex=1/length(numcol)*(factor+1))
	# 	rect(poslegend,ylegend,poslegend+2,ylegend-0.05*maxyplot,col=colr)
		#group of mutation
	}
	ylegend<--maxyplotminus
	text((listxpos[[1]]+listxpos[[length(lensamples)]]+1)/2,ylegend,labcol,cex=1/length(numcol)*factor,pos=1)
# 	rect(poslegend,ylegend,poslegend+1,ylegend-1,col=col)
# 	print(listma)
namesind<-1:length(names)
	for (i in 1:length(lensamples)){
# 		print(listma[[i]])
		#micros
# 		pvect<-0
		text(listxpos[[i]]+0.5,listypos[[i]][4,],namesind[i],cex=1/length(numcol)*(factor-1),pos=3)
		rect(listxpos[[i]],0,listxpos[[i]]+1,listypos[[i]][1,],col=coln[1])
		rect(listxpos[[i]],listypos[[i]][1,],listxpos[[i]]+1,listypos[[i]][2,],col=coln[2])
		rect(listxpos[[i]],listypos[[i]][2,],listxpos[[i]]+1,listypos[[i]][3,],col=coln[3])
		rect(listxpos[[i]],listypos[[i]][3,],listxpos[[i]]+1,listypos[[i]][4,],col=coln[4])
		if (is.null(list$pvaluent)==F){
# 			print("ok")
			pvect<-mapply(getpvaluesNT,listma[[i]][1,],labcol,apply(listma[[i]][9:13,],2,sum),'A',type,i)
# 			print(pvect)
			ypvalue<-listypos[[i]][1,]/2
			xpvalue<-(listxpos[[i]][pvect>0]+listxpos[[i]][pvect>0]+1)/2
			if (length(xpvalue)>0 ){points(xpvalue,ypvalue[pvect>0],pch="*",lwd=.2,col="white")}
			pvect<-mapply(getpvaluesNT,listma[[i]][2,],labcol,apply(listma[[i]][9:13,],2,sum),'T',type,i)
# 			print(pvect)
			xpvalue<-(listxpos[[i]][pvect>0]+listxpos[[i]][pvect>0]+1)/2
			ypvalue<-listypos[[i]][1,]+(listypos[[i]][2,]-listypos[[i]][1,])/2
			if (length(xpvalue)>0){points(xpvalue,ypvalue[pvect>0],pch="*",lwd=.2,col="white")}
			pvect<-mapply(getpvaluesNT,listma[[i]][3,],labcol,apply(listma[[i]][9:13,],2,sum),'C',type,i)
# 			print(pvect)
			xpvalue<-(listxpos[[i]][pvect>0]+listxpos[[i]][pvect>0]+1)/2
			ypvalue<-listypos[[i]][2,]+(listypos[[i]][3,]-listypos[[i]][2,])/2
			if (length(xpvalue)>0){points(xpvalue,ypvalue[pvect>0],pch="*",lwd=.2,col="white")}
			pvect<-mapply(getpvaluesNT,listma[[i]][4,],labcol,apply(listma[[i]][9:13,],2,sum),'G',type,i)
# 			print(pvect)
			xpvalue<-(listxpos[[i]][pvect>0]+listxpos[[i]][pvect>0]+1)/2
			ypvalue<-listypos[[i]][3,]+(listypos[[i]][4,]-listypos[[i]][3,])/2
			if (length(xpvalue)>0){points(xpvalue,ypvalue[pvect>0],pch="*",lwd=.2,col="white")}
			
		}
		#ratios
		colref<-1
		rect(listxpos[[i]],0,listxpos[[i]]+1,-listypos[[i]][iniref,],col=colr[colref])
		
		for (nref in (iniref):(endref-1)){
			colref<-colref+1
			rect(listxpos[[i]],-listypos[[i]][nref,],listxpos[[i]]+1,-listypos[[i]][nref+1,],col=colr[colref])
# 			pvect<-mapply(getpvaluesNT,listma[[i]][1,],labcol,apply(listma[[i]][9:13,],2,sum),'A',type,i)
# 			print(pvect)
# 			ypvalue<-listypos[[i]][1,]/2
# 			xpvalue<-(listxpos[[i]][pvect>0]+listxpos[[i]][pvect>0]+1)/2
# 			if (length(xpvalue)>0 ){points(xpvalue,ypvalue[pvect>0],pch="*",lwd=.2,col="white")}
		}
		
	}
	var<-0.05*maxyplot
	colp<-rainbow(length(lensamples))
	nc<-0
	pcutoff<-0.05
	
	for (m1 in 1:(length(lensamples))){
		nc<-nc+1
		ytopfix<-apply(summaryN,2,max)+nc*var
		
		for (m2 in (m1):length(lensamples)){
					
			
			
			p<-mapply(pvalue,rep(lensamples[m1],ncol(listma[[m1]])),rep(lensamples[m2],ncol(listma[[m1]])),apply(listma[[m1]][9:13,],2,sum),apply(listma[[m2]][9:13,],2,sum))
			if(length(p[p<pcutoff])>0){
# 				ytop<-seq(ytopfix,ytopfix+var*length(p[p<0.05]),var)
				ytop<-ytopfix[p<pcutoff]
				x1<-(listxpos[[m1]][p<pcutoff]+listxpos[[m1]][p<pcutoff]+1)/2
				x2<-(listxpos[[m2]][p<pcutoff]+listxpos[[m2]][p<pcutoff]+1)/2
				y1<-listypos[[m1]][4,p<pcutoff]+var/2
				y2<-listypos[[m2]][4,p<pcutoff]+var/2
# 				ytop<-max(summary)
				segments(x1,y1,x1,ytop,colp[nc])
# 				segments(listxpos[[1]][m1],y1,listxpos[[length(lensamples)]][m1]+1,y1,colp[nc])
				segments(x2,y2,x2,ytop,colp[nc])
# 				segments(listxpos[[1]][m2],y2,listxpos[[length(lensamples)]][m2]+1,y2,colp[nc])
				segments(x1,ytop,x2,ytop,colp[nc])
# 				print(c(m1,m2,p))
				
			}
		}
# 		var<-var+var
	}
	
	dev.off()
	
	doc = newXMLDoc()
	top=newXMLNode("body")
	newXMLNode("h3",attrs=c(align="center"),paste("Studying the variability of ",labelsmut[2]),parent=top)
	div<-newXMLNode("div",attrs=c(style="text-align: center;"), parent = top)
	newXMLNode("img", attrs = c(src = paste(sep="","img.jpg")), parent = div)
	labelsrow<-c("A","U","C","G","A","U","C","G",">80","60-80","60-40","40-20","<20")
	script<-newXMLNode("script", attrs=c(language="JavaScript",src="showtableid.js"),"",parent=top)
	tblbig<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
	trbig<-newXMLNode("tr",  parent = tblbig)
	tdbig<-newXMLNode("td","label" ,parent = trbig)
	tdbig<-newXMLNode("td","name" ,parent = trbig)
	for (ind in namesind){
		trbig<-newXMLNode("tr",  parent = tblbig)
		tdbig<-newXMLNode("td",ind ,parent = trbig)
		tdbig<-newXMLNode("td",names[ind] ,parent = trbig)
	}
	colrnames<-c(">80","60-80","60-40","40-20","<20")
	br<-newXMLNode("p","Percentage of the isomiRs with respect to the corresponding reference miRNAs. ",attrs=c(align="center"),parent=top)
	
	text<-"The variability significance is calculated using the following equation: r=Fv/(Fr+Fv)*100, where Fr is the frequency of the reference sequence and Fv is the frequency of the variant sequence."
	newXMLNode("a","?",attrs=c(title=text,href=""),parent=br)
	newXMLNode("p",parent=top)
	tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = top)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(colall[9:13],function(x) newXMLNode("td", attrs=c(bgcolor=x,style="height:1em;"),parent=tr))
	addChildren(tr, col)
	tr<-newXMLNode("tr",  parent = tbl)
	col=lapply(colrnames,function(x) newXMLNode("td", x,parent=tr))
	addChildren(tr, col)
# 	newXMLNode("p",attrs=c(align="center"),"Abundance of isomiRs (%) with respect
# # to the corresponding reference miRNA",parent=top)
	newXMLNode("p",parent=top)
	
	ptag<-newXMLNode("p",parent = top)
	tblbig<-newXMLNode("table",attrs=c(align= "center",border=0), parent = top)
	trbig<-newXMLNode("tr",  parent = tblbig)
	tdbig<-newXMLNode("td", attrs=c(colspan=length(labcol)),"NTs with respect to
the reference miRNA", parent = trbig)
	trbig<-newXMLNode("tr",  parent = tblbig)
	for (p in numcol){
		
		tdbig<-newXMLNode("td",  parent = trbig)
# 		ptag<-newXMLNode("p",paste("position",p),parent = tdbig)
		tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = tdbig)
		tr<-newXMLNode("tr",  parent = tbl)
		td<-newXMLNode("td",attrs=c(colspan=2),paste("position:",labcol[numcol==p]),parent = tr)	
# 		for (ns in 1:length(names)){
# 			td<-newXMLNode("td",attrs=c(bgcolor="orange"),names[ns], parent = tr)
# 		}
		for (nr in 1:4){
			tr<-newXMLNode("tr",  parent = tbl)
			
# 			for (ns in 1:length(names)){
				index<-nr
				colorfont<-"black"
# 				text<-paste(listmu[[nm]][[index]])
				#text<-paste(sep="","#",p,nr)
                text<-paste(sep="","javascript:loadtable('",p,nr,"')")
				td<-newXMLNode("td",attrs=c(bgcolor=colall[nr],width=2),parent = tr)
				td<-newXMLNode("td",parent = tr)
				a<-newXMLNode("a",attrs=c(href=text),labelsrow[nr],parent=td)
				
# 			}
		}
	
	}
	if (type=="mut"){
	trbig<-newXMLNode("tr",  parent = tblbig)
	tdbig<-newXMLNode("td", attrs=c(colspan=length(labcol)),"NTs in isomiR sequence", parent = trbig)
	trbig<-newXMLNode("tr",  parent = tblbig)
		for (p in numcol){
		
		tdbig<-newXMLNode("td", parent = trbig)
# 		ptag<-newXMLNode("p",paste("position",p),parent = tdbig)
		tbl<-newXMLNode("table",attrs=c(align= "center",border=1), parent = tdbig)
		tr<-newXMLNode("tr",  parent = tbl)
		td<-newXMLNode("td",attrs=c(colspan=2),paste("position:",labcol[numcol==p]),parent = tr)	
# 		for (ns in 1:length(names)){
# 			td<-newXMLNode("td",attrs=c(bgcolor="orange"),names[ns], parent = tr)
# 		}
		for (nr in 1:4){
			tr<-newXMLNode("tr",  parent = tbl)
			
# 			for (ns in 1:length(names)){
				index<-nr
# 				text<-paste(listmu[[nm]][[index]])
				#text<-paste(sep="","#",p,nr+4)
                text<-paste(sep="","javascript:loadtable('",p,nr,"')")
				td<-newXMLNode("td",attrs=c(bgcolor=colall[nr],width=2),parent = tr)
				td<-newXMLNode("td",parent = tr)
				
				a<-newXMLNode("a",attrs=c(href=text),labelsrow[nr],parent=td)
# 			}
		}
		
		}
	}
	div<-newXMLNode("div",attrs=c(id="show",style="text-align: center;"), "test",parent = top)
	for (p in numcol){
		
		newXMLNode("p",parent=top)	
		
		for (nr in numrow){
			iden<-paste(sep="",p,nr)
			
			table2<-matrix(nrow=1,ncol=length(names)+1)
			table2[is.na(table2)]<-0
			table2<-data.frame(table2)
			names(table2)<-c("chr",names)
			table2ratio<-table2
# 			for (nra in 1:5){
			for (ns in 1:length(names)){
				templab<-paste(sep="",nts2[nr])
				index<-txs*ns-txs+indlab[lab==templab]
				temp<-data.frame(chr=(listmu[[p]][[index]]),flag=(listfreqmu[[p]][[index]]))
# 				print(listfreqmu[[p]][[index]])
				tempratio<-data.frame(chr=(listmu[[p]][[index]]),flag=(listmuratio[[p]][[index]]))
				if (nrow(temp)==0){

					temp<-data.frame(chr=0,flag=0)
					names(temp)<-c("chr",names[ns])
					tempratio<-temp
				}
# 				print (nra)
# 				print (temp)
				names(temp)<-c("chr",names[ns])
				names(tempratio)<-c("chr",names[ns])
				if (ns==1){
					table<-temp
					tableratio<-tempratio
				}else{
					table<-merge(temp,table,by="chr",all=TRUE)
					tableratio<-merge(tempratio,tableratio,by="chr",all=TRUE)
				}
				
			}
# 				table<-cbind(rep(nra,nrow(table)),table)
# 				print(table2ratio)
# 				print(tableratio)
# 				names(table)[1]<-"r"
				table2<-rbind(table2,table)
				table2ratio<-rbind(table2ratio,tableratio)
# 			}
# 			print(table)
			table2[is.na(table2)]<-0
			table<-table2[table2$chr!=0,]
			table2ratio[is.na(table2ratio)]<-6
			tableratio<-table2ratio[table2ratio$chr!=0,]
# 			print(table)
			
            tbl<-newXMLNode("table",attrs=c(style="visibility:hidden;",id=iden,align= "center",border=1), parent = top)
			#tbl<-newXMLNode("table",attrs=c(id=iden,align= "center",border=1), parent = top)
			tr<-newXMLNode("tr",parent = tbl)
			td<-newXMLNode("th",attrs=c(colspan=ncol(table)),paste("position:",labcol[numcol==p],"NT:",labelsrow[nr]),parent=tr)
			tr<-newXMLNode("tr",parent = tbl)
			
			td<-newXMLNode("th",attrs=c(bgcolor=colall[nr]),labelsrow[nr],parent=tr)
			col=lapply(names,function(x) newXMLNode("th", x,parent=tr))
			
			addChildren(tr, col)
# 			print("ok")
# 			print(table)
			if (nrow(table)>0){
			
			
			for (nr in 1:nrow(table)){
				tr<-newXMLNode("tr",parent = tbl)
# 				print(table[nr,])
# 				print(table)
# 				print(text)
# 				td<-newXMLNode("td", attrs=c(bgcolor=text,width=2),parent=tr)
				text<-as.character(table[nr,1])
# 				print(table[nr,])
# 				print(tableratio[nr,2:ncol(tableratio)])
				color<-as.character(colall[as.numeric(tableratio[nr,2:ncol(tableratio)]+8)])
# 				print(color)
				color[is.na(color)]<-"white"
				colormi<-"white"
				if (prod(table[nr,2:ncol(table)])>0) {
					colormi<-"#FF6666"
				}
				td<-newXMLNode("td",attrs=c(bgcolor=colormi),text,parent=tr)
				colorfont<-"black"
				if (ncol(table)>2){
					text<-table[nr,2:ncol(table)]
					col=mapply(dotds,text,color,tr)
# 					addChildren(tr, col)
				}else{
					text<-table[nr,2]
					td<-newXMLNode("td",attrs=c(bgcolor=color),text,parent=tr)
				}
				

			}
			}
		}
	
	}
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

if (is.null(list$DB)==F){
	
	options<-paste(sep="",options," `DB` LIKE '",list$DB,"' AND")
}
if (is.null(list$locus)==F){
	
	options<-paste(sep="",options," `chr` LIKE '%",list$locus,"%' AND")
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
typemut<-0
labelsmut<-0
iniref<-9
endref<-13
type<-0
listmu<-vector("list",length=25)
listfreqmu<-vector("list",length=25)
listmuratio<-vector("list",length=25)

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
	labelsmut<-c(labelsmut,"nt substitution")
	numcol<-c(ini:size)
	labcol<-c(ini:size)
	pospvalue<-c(ini:size)
	iniref<-5
	endref<-8
	type<-"mut"
}
nts<-c("A","T","C","G")
ntspos<-c(1,2,3,4)
nts2<-c("A","T","C","G","mA","mT","mC","mG")
lab<-0
for (nt in nts2){

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
}

ns<-1
freqNTHairpin<-vector("list",length=length(list$group1))
listsamples<-vector("list",length=length(list$group1))
normalize<-vector("list",length=length(list$group1))
lensamples<-0
listprofile<-vector("list",length=length(list$group1))
listperfect<-vector("list",length=length(list$group1))
listma<-vector("list",length=length(list$group1))
ma<-matrix(ncol=25,nrow=13)
ma[is.na(ma)]<-0


count<-0
checktrimming5<-0
checktrimming3<-0
for (s in list$group1){
	s<-paste(sep="",list$project,"`.`",s)
	if (type=="mut"){
		
		query<-paste(sep="","select`chr`,`seq`,`mut`,`start`,`end`,`trimmed5` from `",s,"` ",options," AND `amb`=1 AND (`trimmed5` like '0' OR  `trimmed5` like 'na') AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+' AND `mut` NOT REGEXP '^[0-9]+[ACTG]+[0-9]+' GROUP BY `chr`;")
		
		rs <- dbSendQuery(con,query) 
		temp <- as.data.frame(fetch(rs))
# 		print((temp))
		table<-matrix(nrow=nrow(temp),ncol=size)
		table[,]<-mapply(dorefseq,temp$seq,temp$mut,temp$trimmed5,temp$end-temp$start+1,ini,size)
# 		print ("ok")
		freqNTHairpin[[ns]]<-table
# 		print(freqNTHairpin[[ns]])
		indlenreal<-1:length(freqNTHairpin[[ns]][1,])
		indtrlen<-1:length(indlenreal)
		
	}else if (type!="3 addition"){
	
		query<-paste(sep="","select`chr`,`seq`,`trimmed5ref`,`trimmed3ref` from `",s,"` ",options," AND `amb`=1 AND (`trimmed5` like '0' OR  `trimmed5` like 'na') GROUP BY `chr`;")
		
		rs <- dbSendQuery(con,query) 
		temp <- as.data.frame(fetch(rs))
		if (length(temp$chr[temp$trimmed5ref!="na"])==0 & nrow(temp)>0){
			checktrimming5<-1
		}
		if (length(temp$chr[temp$trimmed3ref!="na"])==0 & nrow(temp)>0){
			checktrimming3<-1
		}
		table<-matrix(nrow=nrow(temp),ncol=size-ini+1)
		if (type=="trimmed5"){table[,]<-mapply(dorefseqtrm,temp$trimmed5ref,size)}
		if (type=="trimmed3"){table[,]<-mapply(dorefseqtrm,temp$trimmed3ref,size)}
		freqNTHairpin[[ns]]<-table
# 		print(freqNTHairpin[[ns]][1:10,])
		lenreal<-length(freqNTHairpin[[ns]][1,])
		indlenreal<-c(-(lenreal/2):-1,1:(lenreal/2))
		indtrlen<-1:length(indlenreal)
	}
	
	query<-paste(sep="","select `id`,`chr`,`freq`,`trimmed5`,`trimmed3`,`addition3`,`mut`,`len` from `",s,"` ",options," AND `amb`=1 AND  `mut` NOT REGEXP '^0[ACTG]+' AND  `mut` NOT REGEXP '^-[0-9]+[ACTG]+' AND `mut` NOT REGEXP '^[0-9]+[ACTG]+[0-9]+';")
	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	listsamples[[ns]]<-temp
	lensamples[ns]<-length(unique(listsamples[[ns]][,2]))
# 	if (is.null(list$premapped)==T){
		listperfect[[ns]]<-temp[temp$trimmed5=="0" & temp$trimmed5!="na" & temp$trimmed3=="0" & temp$trimmed3!="na" & temp$addition3=="0" & temp$addition3!="na" & temp$mut=="0" & temp$mut!="na",2:3]
	normalize[[ns]]<-1:25
	for (np in numcol){
# 		normalize[[ns]][np]<-length(unique(temp$chr[temp$trimmed5=="0" & temp$trimmed5!="na" & temp$trimmed3=="0" & temp$trimmed3!="na" & temp$addition3=="0" & temp$addition3!="na" & temp$mut=="0" & temp$mut!="na" & temp$len>=np]))
		normalize[[ns]][np]<-length(unique(temp$chr[temp$len>=np]))
	}
# 	print (normalize)
	if (ns==1){
		listmi<-unique(listsamples[[ns]][,2])
	}else{
		listmi<-unique(c(listmi,listsamples[[ns]][,2]))
	}
	listprofile[[ns]]<-0
	###########################WARNINGGGGGGGGGGG

	listma[[ns]]<-ma
	ns<-ns+1
}
# q()
if (is.null(list$lmi)==F){
	query<-paste(sep="","select `id`,`len` from `DB",list$DB,"`;")
	print(query)
	rs <- dbSendQuery(con,query) 
	lenmi <- as.data.frame(fetch(rs))
	listmi<-intersect(listmi,lenmi$id[lenmi$len==list$lmi])

}

# ma<-as.data.frame(ma)
# print(ma)

tempmaind<-1:25
allisomirs<-data.frame(id=0,freq=0)
# total<-0
if (length(listmi)>0 ){

for (chr in listmi){
# 	print (chr)
	for (ns in 1:length(list$group1)){
# 		print (ns)
# 		if (type!="addition3"){
# 			listprofile[[ns]]<-append(listprofile[[ns]],indfreqNT[freqNTHairpin[,1]==chr])
# 		}
		tempperfect<-listperfect[[ns]]
# 		print (tempperfect)
		if (length(tempperfect$freq[tempperfect$chr==chr])==0){
			freqperfect<-0
		}else{
			freqperfect<-sum(tempperfect$freq[tempperfect$chr==chr])
		}
		for (nm in typemut[2:length(typemut)]){
# 			print(nm)
			tempdata<-listsamples[[ns]]
			if (length(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])>0){
				
				tempvalue<-sum(tempdata$freq[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na"])
				minvalue<-min(freqperfect,tempvalue,na.rm=T)
# 				print(c(freqperfect,tempvalue))
				temp<-tempdata[tempdata$chr==chr & tempdata[,nm]!="0" & tempdata[,nm]!="na",]
# 				if (chr=="hsa-miR-192"){
# 					print(c("before",nrow(temp)))
# 					print(temp[1:10,])
# 				}
				temp<-filter(temp,minvalue*cof,nm,"filter",freqperfect,list$error)
# 				if (chr=="hsa-miR-192"){
# 					print(nrow(temp))
# 				}
				if (length(temp$id)>0){
				freqperfect<-freqperfect+(tempvalue-sum(temp$freq,na.rm=T))
				tempma<-matrix(nrow=8,ncol=25)
				tempma[is.na(tempma)]<-0
				for (i in 1:length(temp$id)){
# 				tempnisomirs<-temp[,c(1,3)]
# 				allisomirs<-merge(temp,allisomirs,by=1,all=T)
# 				allisomirs<-allisomirs[,1:2]		
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
# 						if (chr=="hsa-miR-140-3p"){
# 							print(temp$freq[i])
# 						}
					}else if (nm==5){
						mutation<-unlist(strsplit(temp[i,nm],""))
# 						print(mutation)
						nt<-mutation[length(mutation)]
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
						mutation<-unlist(strsplit(gsub("[0-9]+","",temp[i,nm]),""))
						pos<-as.numeric(strsplit(gsub("[ATGC]+"," ",temp[i,nm])," "))
# 						if(pos<)
# 						print(mutation)
						countnt<--1
						for (posmutmir in pos){
							countnt<-countnt+2
							nt1<-mutation[countnt]
							nt2<-mutation[countnt+1]
							tempma[ntspos[nts==nt2],posmutmir]<-tempma[ntspos[nts==nt2],posmutmir]+temp$freq[i]
							tempma[ntspos[nts==nt1]+4,posmutmir]<-tempma[ntspos[nts==nt1]+4,posmutmir]+temp$freq[i]
						}
# 						total<-total+nrow(temp)
					}
					
				}
# 				print(chr)
# 				print(tempma)
# 				if (chr=="hsa-miR-192"){
# 					print(tempma[,numcol])
# 				}
				ma<-listma[[ns]]
				maxvalue<-apply(tempma,2,max)
				for (maxn in numcol){
					ratio<-freqperfect/(freqperfect+maxvalue[maxn])
					if (is.na(ratio)==TRUE){
						ratio<-1
					}
					if (is.na(ratio)==FALSE & maxvalue[maxn]>0){
						value<-as.numeric(cut(ratio,breaks=c(-1,seq(0.2,1,0.2)),labels=1:5,right=T))
						ma[value+8,maxn]<-ma[value+8,maxn]+1
# 						index<-13*ns-13+value+8
# 						print 
# 						listmu[[maxn]][[index]]<-c(listmu[[maxn]][[index]],chr)
# 						listfreqmu[[maxn]][[index]]<-c(listfreqmu[[maxn]][[index]],maxvalue[maxn])
# 						print (c(value,ratio))
					}
				}
				
				
				
				for (typent in 1:8){
					if (length(tempma[tempma[typent,numcol]>0])>0){
					ma[typent,tempmaind[tempma[typent,]>0]]<-ma[typent,tempmaind[tempma[typent,]>0]]+1
# 					print(tempma[1:8,10:15])
# 					print(ma[1:8,10:15])
					}
					
					for (p in numcol){
						
						if (tempma[typent,p]>0){
						count<-count+1
						ratio<-freqperfect/(freqperfect+tempma[typent,p])
						value<-as.numeric(cut(ratio,breaks=c(-1,seq(0.2,1,0.2)),labels=1:5,right=T))
# 						
# 						pre<-""
# 						if(typent>4){
# 							pre<-"m"
# 						}
						templab<-paste(sep="",nts2[typent])
						index<-txs*ns-txs+indlab[lab==templab]
# 						print(c(templab))
						listmu[[p]][[index]]<-c(listmu[[p]][[index]],chr)
						listfreqmu[[p]][[index]]<-c(listfreqmu[[p]][[index]],tempma[typent,p])
						listmuratio[[p]][[index]]<-c(listmuratio[[p]][[index]],value)
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

}#for
}#if

# total
#pospvalue
# q()
# print(count)
writexml(listma,lensamples,list$group1,labelsmut,labcol,numcol,pospvalue,type,iniref,endref)

dbDisconnect(con) 
