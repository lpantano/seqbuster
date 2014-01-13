
chosecolor<-function(n){

	if (n<10){
		col<-c("#00FFCC","#00CCCC","#0099CC","#0066CC","#0033CC","#0000CC","#000099","#000066","#000033","#000000")
	}else{
		col=rainbow(n)
	}
	if (n==1){
		col<-col[1]
	}else{
		col<-col[1:n]
	}
	return (col)

}

pvalue<-function(sum1,sum2,obs1,obs2) {
if (is.na(obs1)){obs1=0}
if (is.na(obs2)){obs2=0}
if (sum(obs1,obs2)!=0){
	m<-matrix(ncol=2,nrow=2)
	m[1,]<-c(obs1,obs2)
	m[2,]<-c(sum1,sum2)
	f<-chisq.test(m)
	p<-f$p.value
}else{

	p<-1
}


return(p)


}

BHcorrection<-function(p,count,rank){

	q<-p*count/rank
	return (q)
}




ztest<-function(x1,x2,n1,n2){

	ps=(x1+x2)/(n1+n2)
	qs=1-ps
	ns=(1/n1+1/n2)
	p1=x1/n1
	p2=x2/n2
	num=p1-p2
	div=sqrt(ps*qs*ns)
	delta=abs(num/div)
	p<-2*dnorm(delta)
	return(p)
}


scalefreq<-function(freqs,sc){
# 	print (sum(freqs))
	newfreqs<-round(freqs/sum(freqs)*sc)
	return (newfreqs)
}


logfreq<-function(freqs,l){
	if (l=="ln"){
		l=exp(1)
	}else if(l=="log10"){
		l=10
	}else if(l=="log2"){
		l=2
	}
	
# 	print(freqs)
	newfreqs<-(log(freqs,base=l))
	newfreqs[is.infinite(newfreqs)]<-0
# 	print(newfreqs)
	return (newfreqs)
}

discartfreq<-function(freqs,qmax,qmin){
# 	print(freqs)
	q<-as.numeric(quantile(as.numeric(unlist(unique(freqs[,2]))),probs=c(qmin/100,qmax/100)))
# 	print(q[1])
	newfreqs<-freqs[freqs[,2]>=q[1] & freqs[,2]<=q[2],]
# 	print(q)
	return (newfreqs)
}

getquantilefreq<-function(freqs,qmax,qmin){

	q<-as.numeric(quantile(as.numeric(unlist(unique(freqs[,2]))),probs=c(qmin/100,qmax/100)))
	newfreqs<-freqs[freqs[,2]>=q[1] & freqs[,2]<=q[2],]
# 	print(q)
# 	print(newfreqs)
	return (newfreqs)
}

centernorm<-function(freqs,op){
	newfreqs<-freqs
	if (op=="mean"){
		newfreqs[,2]<-freqs[,2]/mean(freqs[,2])
	}else if (op=="median"){
		newfreqs[,2]<-freqs[,2]/median(freqs[,2])
	}else if (op=="moda"){
		newfreqs[,2]<-freqs[,2]/moda(freqs[,2])
	}else if (op=="max"){
		newfreqs[,2]<-freqs[,2]/max(freqs[,2])
	}else if (op=="min"){
		sortfreq<-sort(freqs[,2])
# 		print((sortfreq[1:10]))
		newfreqs[,2]<-freqs[,2]/median(sortfreq[1:10])
	}else if (op=="quantile"){
		q<-as.numeric(quantile(unlist(unique(freqs[,2])),probs=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)))
		q[1]<-0
		for (qn in 2:length(q)){
			m<-mean(freqs[freqs[,2]>q[qn-1] & freqs[,2]<=q[qn],2])
			newfreqs[newfreqs[,2]>q[qn-1] & newfreqs[,2]<=q[qn],2]<-freqs[freqs[,2]>q[qn-1] & freqs[,2]<=q[qn],2]/m
		}
	}else if (op=="houskeeping"){
		freqn<-freqs[freqs[,1]==houskeeping,2]
		newfreqs[,2]<-freqs[,2]/freqn
	}
	newfreqs[,2]<-round(newfreqs[,2])
	return (newfreqs)
}



normbyhyper<-function(freqs,pop,qmax,qmin){

	newfreqs<-discartfreq(freqs,qmax,qmin)
	tend<-sum(newfreqs[,2])
	
	for (i in 1:length(newfreqs[,2])){

		p<-newfreqs[i,2]/tend
		m<-p*pop
		n<-pop-p*pop

		v<-qhyper(p,m,n,pop)
		newfreqs[i,2]<-v
	}
	return(newfreqs)

}


checkchr<-function(chr){
	r<-0
	if (length(list$loci[list$loci==chr])>0){
		r<-1
	}
	return (r)

}
getcommon<-function(data,n){
	scale<-round(log(data$freq,base=2))
	freq<-scale[n]
	data<-data[log(data$freq,base=2)>=freq,]
	return (data)
}

applynorm<-function(tempn,scaleval,q2val,q1val,typetrval,typetdval){
#  	print(c(scaleval,q1val,q2val,typetrval,typetdval))
	
	if (as.numeric(q1val)<100 | as.numeric(q2val)>0 ){
		tempn<-discartfreq(tempn,as.numeric(q1val),as.numeric(q2val))
	}
	
	if(typetdval!=""){
		tempn<-centernorm(tempn,typetdval)
	}
	if (scaleval!="" ){
		tempn$freq<-scalefreq(tempn$freq,as.numeric(scaleval))
	}
	if (typetrval!=""){
		tempn$freq<-logfreq(tempn$freq,(typetrval))
	}
	

	return(tempn)
}
library(RMySQL)

ttestfun<-function(gd1,gd2){
	vec<-unlist(as.vector(append(as.vector(gd1),as.vector(gd2))))
	n1<-length(gd1)
	n2<-length(gd2)
	p<-t.test(gd1,gd2)
	pval<-p$p.value
	sigma<-sqrt((var(gd1)+var(gd2))/2)
	d<-abs(mean(gd1)-mean(gd2))/sigma
	es<-d*sqrt((n1*n2)/(n1+n2))
	return (c(pval,es))
}

anovafun<-function(gd1,gd2){
	vec<-unlist(as.vector(append(as.vector(gd1),as.vector(gd2))))
	podwt <- data.frame(wt=vec,treat=c(rep("G1",length(gd1)),rep("G2",length(gd2))))
	fitpodwt <- lm(wt~treat, data=podwt)
	p<-anova(fitpodwt)
	pval<-unlist(p[5])[1]
	ss<-p[,2][2]
	p1<-length(gd1)/sum(length(gd1),length(gd2))
	p2<-length(gd2)/sum(length(gd1),length(gd2))
	m1<-mean(gd1)
	m2<-mean(gd2)
	gm<-sum(gd1,gd2)/sum(length(gd1),length(gd2))
	es<-sqrt(((p1*(m1-gm)*(m1-gm))+(p2*(m2-gm)*(m2-gm)))/ss)
	return (c(pval,es))
}

difexpseq<-function(stat,scaleval,q1val,q2val,samplesnames1,samplesnames2,projectval,typetrval,typetval,localhost,user,pssw,port){
MySQL(max.con = 16, fetch.default.rec = 5000000, force.reload = TRUE)

m <- dbDriver("MySQL")



if (port==0){
	con <- dbConnect(m,host=localhost,user=user,db=projectval,password=pssw)
}else{
	con <- dbConnect(m,host=localhost,user=user,db=projectval,password=pssw,port=port)
}

cs<-0.5
cb<-1.5
pval<-0.05

infogroup<-vector()

scale<-as.numeric(scaleval)
ns<-0
listsamples1<-vector("list",length=length(samplesnames1))
max1<-1:length(samplesnames1)
table<-data.frame(id=0,freq=0)
listseqnames<-data.frame(seq='0',name='0')
for (s in samplesnames1){
	ns<-ns+1
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `idu`,`freq` from `",projectval,"`.`",s,"clusmap` where`idu`>0 group by `idu`;")
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",projectval,"`.`",s,"clusraw` where `id`>0;")
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
		table<-applynorm(temp,scaleval,q1val,q2val,typetrval,"")
		if (nrow(table)==0){
			table<-data.frame(idu=0,freq=0)
			names(table)<-c('idu',s)
		}
    if (ns==1){
		all<-table
		names(all)<-c('idu',s)
	
	}else{
		table<-table
		names(table)<-c('idu',s)
		all<-merge(all,table,by="idu",all=TRUE)
		
	}
}
#calculate pvalue intra groups
all[is.na(all)]<-0
table1<-all
ncol1<-ncol(table1)
ns<-0
listsamples2<-vector("list",length=length(samplesnames2))
max2<-1:length(samplesnames2)
#Load samples inside groups
for (s in samplesnames2){
	ns<-ns+1
	infogroup<-append(infogroup,paste(sep="","G2:",s))

	query<-paste(sep="","select `idu`,`freq` from `",projectval,"`.`",s,"clusmap` where `idu`>0 group by `idu`;")
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",projectval,"`.`",s,"clusraw` where `id`>0;")
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
		table<-applynorm(temp,scaleval,q1val,q2val,typetrval,"")
	if (nrow(table)==0){
			table<-data.frame(idu=0,freq=0)
			names(table)<-c('idu',s)
		}
	if (ns==1){
		all<-table
		names(all)<-c('idu',s)
		
	}else{
		table<-table
		names(table)<-c('idu',s)
		all<-merge(all,table,by="idu",all=TRUE)
		
	}
}

#calculate pvalue intra groups
all[is.na(all)]<-0
table2<-all
ncol2<-ncol(table2)
nstart2<-ncol1+1
###########join tables

table<-merge(table1,table2,by=1,all=TRUE)
table[is.na(table)]<-0


for (i in 1:nrow(table)){

	if (stat=="Anova"){

		t<-anovafun(table[i,2:ncol1],table[i,nstart2:ncol(table)])
		table$pval[i]<-t[1]
		table$ef<-t[2]
	}else{
		t<-ttestfun(table[i,2:ncol1],table[i,nstart2:ncol(table)])
		table$pval[i]<-t[1]
		table$ef<-t[2]
	}
}
sort<-sort(table$pval,index.return=T)
ind<-1:nrow(table)
order<-unlist(lapply(ind,function (x) ind[sort$ix==x]))
table$qvalue<-mapply(BHcorrection,table$pval,order,nrow(table))



 query<-paste(sep="","DROP TABLE IF EXISTS `",projectval,"`.`difexp`;")
 rs <- dbSendQuery(con,query) 

createtable<-paste(sep="","CREATE TABLE `",projectval,"`.`difexp` (`id` int unsigned, ")
for (c in 2:(ncol(table)-3)){
	createtable<-paste(sep="",createtable,"`",names(table)[c],"` DOUBLE(7,2), ")
}
createtable<-paste(sep="",createtable," `pval`  DOUBLE(8,2),qval DOUBLE(8,2),effect  DOUBLE(8,2),INDEX (id));")
 rs <- dbSendQuery(con,createtable) 

for (r in 1:nrow(table)){
 	
	
	loadtable<-paste(sep="","INSERT INTO `",projectval,"`.`difexp` VALUES (")
	for (c in 1:(ncol(table)-1)){
			loadtable<-paste(sep="",loadtable,table[r,c],", ")
	}
	loadtable<-paste(sep="",loadtable,table[r,ncol(table)]," );")

 	rs <- dbSendQuery(con,loadtable) 
	
}

d<-getwd()

write.table(table,paste(sep="",d,"/","temp/df.done"))
dbDisconnect(con) 

}
d<-getwd()
name<-paste(sep="",d,"/","temp/args.txt")
d<-read.table(name)

opt<-as.vector(d[,1])
# opt
g1<-unlist(strsplit(opt[2],":"))
# g1
g2<-unlist(strsplit(opt[3],":"))
# g2
source("Rscripts/R/db.R")
 difexpseq(opt[1],1000000,0,100,g1,g2,opt[4],"log2","scale",hostname,username,pssw,port)
# difexpseq("Ztest","BH","1000000",0,100,"test2","test2","testpro","","scale","localhost","lpantano","sqllorena")

