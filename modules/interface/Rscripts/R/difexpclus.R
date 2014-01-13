
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

binomial<-function(n1,n2,m1,m2){

	n<-min(n1,n2)
	nm<-max(n1,n2)
	if (n1<50| n2<50){
		p<-sum(dbinom(1:n,n1+n2,0.5))
# 		p<-10
	}else{
		total<-n1+n2
		p<-sum(dnorm(1:n,mean=total*0.5,sd=sqrt(total*0.25)))
	}
	return(p)
}

statistic<-function(n1,n2,ratio){
	if(n1==0){n1=1}
	if(n2==0){n2=1}
	if(n1==n2){
		p<-1
	}else if (n1+n2<170){
		p<-bayesian2(n1,n2,ratio)
	}else if(n1+n2>=170){
		p<-bayesian(n1,n2,ratio)
	}

return (p)
}


bayesian<-function(n1,n2,ratio){
ny<-max(n1,n2)
nx<-min(n1,n2)
# print (c(n1,n2))
x<-seq(1,nx,1)
xy<-seq(ny+1,ny+nx,1)
r<-ratio
r2<-r+1
vr<-rep(r,nx)
vr2<-rep(r,ny-nx)
divr<-rep(r2,nx)
divr2<-rep(r2,ny-nx)
divr3<-1/rep(r2,ny+1-(ny-nx))
divn<-xy/x

comp1<-vr*divn/(divr)

comp2<-vr2/divr2

nmax<-max(nx,ny-nx,ny+1-(ny-nx))
total<-1

mins<-min(nx,ny-nx)
simp<-comp1[1:mins]*comp2[1:mins]
if (mins==nx & length(comp2)>mins){
	simp<-append(simp,comp2[(mins+1):length(comp2)])
}else if (length(comp1)>mins){
	simp<-append(simp,comp1[(mins+1):length(comp1)])
}
mins<-min(length(simp),ny+1-(ny-nx))
simp2<-simp[1:mins]*divr3[1:mins]
if (mins==length(simp) & mins!=length(divr3)){
	simp2<-append(simp2,divr3[(mins+1):length(divr3)])
}else if (mins!=length(simp) & mins==length(divr3)){
	simp2<-append(simp2,simp[(mins+1):length(simp)])
}
temp<-simp2
c<-0
while (length(temp)>1){
c<-c+1
	if (round(length(temp)/2)!=length(temp)/2){
		temp[1]<-temp[1]*temp[length(temp)]
		temp<-temp[1:(length(temp)-1)]
	}
		left<-length(temp)/2
		end<-length(temp)
	
	temp1<-temp[1:left]*temp[(left+1):end]
	temp<-temp1
}



total<-temp
return (total)


}


bayesian2<-function(n1,n2,ratio){
	num<-factorial(n1+n2)
	expn<-n1+n2+1
	den<-factorial(n1)*factorial(n2)*((1+ratio)^expn)
	p<-(ratio^n2)*num/den
	if (num/den==0){
		p<-0
	}
	return (p)
}

filter<-function(data,mv,nmv,type,freq,error){

temp<-data
if (type=="filter"){
# 	print(nmv)
# 	print(temp)
	temp<-data[data$freq>=mv,]
# 	print(temp)
}

if (nmv==7 & error=="zscore" & nrow(temp)>0){
	tempn<-temp[temp$mut!="0",]
	tempnm<-temp[temp$mut=="0",]
	if (nrow(tempn)>0){
		tempn<-zscore(tempn,freq,nmv)
	}
	temp<-rbind(tempn,tempnm)
}

# print(nrow(temp))
# print((temp))
return(temp)

}

bayeschange<-function(data,freq,nmv){
pr<-1/15000/4+0.0005
ppos<-c(0.005,rep(0.003,13),0.004,rep(0.005,3),rep(0.0055,2),rep(0.01,3),rep(0.015,2),rep(0.02,2),rep(0.015,2),0.025,0.03,0.035,0.0375)
pntraw<-matrix(ncol=4,nrow=4)
#cols=from rows=into
pntraw[1,]<-c(0,0.13,0.04,0.08)
pntraw[2,]<-c(0.14,0,0.08,0.15)
pntraw[3,]<-c(0.05,0.02,0,0.09)
pntraw[4,]<-c(0.05,0.04,0.12,0)
#norm<-
pnt<-pntraw/colSums(pntraw)

nt<-c("A","T","C","G")
ntind<-1:4

filter<-data[1,]
# print(data)
m<-ncol(data)
# print ("inside")
# print (nmv)
FDR<-(0.05*(m+1))/(2*m)
# print (FDR)
for (i in 1:nrow(data)){
	mutation<-unlist(strsplit(gsub("[0-9]+","",data[i,nmv]),""))
	ntinto<-mutation[1]
	ntfrom<-mutation[2]
	pos<-as.numeric(gsub("[ATGC]+","",data[i,nmv]))
	pH<-pnt[ntind[ntinto==nt],ntind[ntfrom==nt]]*ppos[pos]
# 	print (pr)
	freqmut<-data$freq[i]
	total<-sum(data$freq,freq)
	pDATA<-freqmut/total
	pDATAH<-pDATA*pH
	if (pH>1){
		p<-0
	}else{
		num<-pDATAH
		div<-pDATA
		r<-num/div
		
		p<-r
		
	}
# 	print (p)
	
	if (p<=FDR){
		filter<-rbind(filter,data[i,])

	}else{
# 		print(freqmut)
	}
}
# print(filter)
if (nrow(filter)==1){
	filter<-filter[filter$id==0,]
}else{
	filter<-filter[2:nrow(filter),]
}
# print(filter)

return(filter)


}

ztest<-function(x1,x2,n1,n2){
	p<-1
	if (x1>0 & x2>0){
		ps=(x1+x2)/(n1+n2)
		qs=1-ps
		ns=(1/n1+1/n2)
		p1=x1/n1
		p2=x2/n2
		num=p1-p2
		div=sqrt(ps*qs*ns)
		delta=abs(num/div)
		p<-2*dnorm(delta)
	}
	return(p)
}

normtest<-function(x1,x2,n1,n2){
	
	p<-pnorm(x1,mean=x2,sd=x2/3)
	if(p>0.9){
		p<-1-p
	}
	return(p)
}



fisher<-function(x1,x2,n1,n2){
	ma<-matrix(ncol=2,nrow=2)
	ma[1,]<-c(x1,x2)
	ma[2,]<-c(n1-x1,n2-x2)
	t<-fisher.test(ma)
	p<-t$p.value
	return(p)
}

zscore<-function(data,freq,nmv){
pr<-1/15000/4+0.0005
ppos<-c(0.005,rep(0.003,13),0.004,rep(0.005,3),rep(0.0055,2),rep(0.01,3),rep(0.015,2),rep(0.02,2),rep(0.015,2),0.025,0.03,0.035,0.0375)
pntraw<-matrix(ncol=4,nrow=4)
#cols=from rows=into
pntraw[1,]<-c(0,0.13,0.04,0.08)
pntraw[2,]<-c(0.14,0,0.08,0.15)
pntraw[3,]<-c(0.05,0.02,0,0.09)
pntraw[4,]<-c(0.05,0.04,0.12,0)
#norm<-
pnt<-pntraw/colSums(pntraw)

nt<-c("A","T","C","G")
ntind<-1:4

filter<-data[1,]
# print(data)
m<-ncol(data)
# print ("inside")
# print (nmv)
FDR<-(0.05*(m+1))/(2*m)
# print (FDR)
for (i in 1:nrow(data)){
#  	print (data[i,])
	mutation<-unlist(strsplit(gsub("[0-9]+","",data[i,nmv]),""))
	ntinto<-mutation[1]
	ntfrom<-mutation[2]
	pos<-as.numeric(gsub("[ATGC]+","",data[i,nmv]))
	pr<-pnt[ntind[ntinto==nt],ntind[ntfrom==nt]]*ppos[pos]
# 	print (pr)
	freqmut<-data$freq[i]
	total<-sum(data$freq,freq)
	pt<-freqmut/total
	
	c1<-pr*total
	c2<-(1-pr)*total
	if (c1<5 | c2<5){
		p<-sum(dbinom(freqmut:total,total,pr))
		
	}else{
		
			num<-pt-pr
			div<-sqrt(pt*(1-pt)/total)
			r<-num/div
			
			p<-1-pnorm(r)
			
		
# 	
	}
	
# 	print (p)
	if (p<=FDR){
		filter<-rbind(filter,data[i,])

	}else{
# 		print(freqmut)
	}
	
}
# print(filter)
if (nrow(filter)==1){
	filter<-filter[filter$id==0,]
}else{
	filter<-filter[2:nrow(filter),]
}
# print(filter)

return(filter)

}




scalefreq<-function(freqs,sc){
#	print (sum(freqs))
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
# 	print(newfreqs)
	newfreqs[,2]<-round(newfreqs[,2])
	return (newfreqs)
}



normbyhyper<-function(freqs,pop,qmax,qmin){

	#v<-dbinom(,1000000,p)*k
# 	uf<-unlist(unique(freqs[,2]))
# 	print(freqs[1:10,])
	newfreqs<-discartfreq(freqs,qmax,qmin)
	#newfreqs2<-getquantilefreq(freqs,qmax,qmin)
# 	newfreqs<-freqs[freqs[,2]<0.2*sum(freqs[,2]),]
# 	q<-as.numeric(quantile(uf,probs=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)))
	tend<-sum(newfreqs[,2])
	
	for (i in 1:length(newfreqs[,2])){

		p<-newfreqs[i,2]/tend
		m<-p*pop
		n<-pop-p*pop
		#k<-
		v<-qhyper(p,m,n,pop)
		newfreqs[i,2]<-v
	}
# 	print(min(newfreqs[,2]))
	#newfreqs<-rbind(newfreqs,newfreqs2)
	return(newfreqs)

}


checkchr<-function(chr){
	r<-0
# 	print ("kk")
	if (length(list$loci[list$loci==chr])>0){
		r<-1
	}
	return (r)

}
getcommon<-function(data,n){
# 	print(n)
	scale<-round(log(data$freq,base=2))
	freq<-scale[n]
#	print(freq)
# 	print(nrow(data))
	data<-data[log(data$freq,base=2)>=freq,]
# 	print(nrow(data))
	return (data)
}
cs<-0.4
cb<-1.5
pval<-0.05

numratiofunc<-function(v){

	n<-length(v[v<cb & v>=cs])
	return(n)
}


numpvalfunc<-function(v){
	
	n<-length(v[v<=pval])
	return(n)
}

numexpfunc<-function(v){
	
	n<-length(v[v<=1])
	return(n)
}


applynorm<-function(tempn,scaleval,q2val,q1val,typetrval,typetdval){
#  	print(c(scaleval,q1val,q2val,typetrval,typetdval))
	
	if (as.numeric(q1val)<100 | as.numeric(q2val)>0 ){
		tempn<-discartfreq(tempn,as.numeric(q1val),as.numeric(q2val))
# 		q<-quantile(tempn[,2],probs=c(10/100,90/100))
# 		print(tempn[1:10,])
		print("quantile")
# 		print(tempn[1:10,])
	}
	
	if(typetdval!=""){
		print("center")
#  		print(tempn[1:10,])
		tempn<-centernorm(tempn,typetdval)
#  		print(tempn[1:10,])
	}
	if (scaleval!="" ){
		print("scale")
# 		print(tempn[1:10,])
		tempn$freq<-scalefreq(tempn$freq,as.numeric(scaleval))
# 		print(tempn[1:10,])
	}
	if (typetrval!=""){
		print("trans")
		tempn$freq<-logfreq(tempn$freq,(typetrval))
# 		print(tempn[1:10,])
	}
	
	
# 	if (list$norm=="hyper"){
# 		tempn<-normbyhyper(tempn,as.numeric(list$scale),as.numeric(list$qmax),as.numeric(list$qmin))
# # 		print(tempn[1:10,])
# 	}
	
	
# 	print(tempn[1:10,])
	
	
	
	return(tempn)
}
library(RMySQL)

# Ztest,BH,1000000,0,100,test2:,test2:,testpro,scale,localhost,lpantano,sqllorena

difexpseq<-function(typemetval,corval,scaleval,q1val,q2val,samplesnames1,samplesnames2,projectval,typetrval,typetdval,localhost,user,pssw,port,outl){

# print (c(typemetval,corval,scaleval,q1val,q2val,samplesnames1,samplesnames2,projectval,typetrval,typetdval,localhost,user))

MySQL(max.con = 16, fetch.default.rec = 5000000, force.reload = TRUE)

m <- dbDriver("MySQL")


	
if (port==0){
	con <- dbConnect(m,host=localhost,user=user,db=projectval,password=pssw)
}else{
	con <- dbConnect(m,host=localhost,user=user,db=projectval,password=pssw,port=port)
}

	query<-paste(sep="","DROP TABLE IF EXISTS `difexp`;")
	rs <- dbSendQuery(con,query) 

cs<-0.5
# cs<-as.numeric(list$cs)
cb<-1.5
# cb<-as.numeric(list$cb)
pval<-0.05
# pval<-as.numeric(list$pval)
outl<-as.numeric(outl)	
	
infogroup<-vector()

scale<-as.numeric(scaleval)
ns<-0
listsamples1<-vector("list",length=length(samplesnames1))
max1<-1:length(samplesnames1)
table<-data.frame(id=0,freq=0)
listseqnames<-data.frame(seq='0',name='0')
#Load samples inside groups
for (s in samplesnames1){
	#s<-paste(sep="",s,projectval)
	ns<-ns+1
	print(s)
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `idu`,`freq` from `",projectval,"`.`",s,"clusmap` where`idu`>0 group by `idu`;")
#  	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",projectval,"`.`",s,"clusraw` where `id`>0;")
# 	print(query)
	rs <- dbSendQuery(con,query) 
	max1[ns]<- unlist(as.vector(fetch(rs)))
# 	if (is.null(list$ref)==F){
# 		temp<-getcommon(temp)
# 	}
  		print(temp[1:5,])
		table<-applynorm(temp[,c(1,2)],scaleval,q1val,q2val,typetrval,typetdval)
#  		print(temp[1:10,c(1,2)])
		table<-merge(temp,table,by=1,all=FALSE)
# 		print(table[1:10,c(1,3)])
	if (ns==1){
		all<-table[,c(1,3)]
		names(all)<-c('id',samplesnames1[ns])
# 		print(all[1:10,])
		#tableall<-temp
		
	}else{
		table<-table[,c(1,3)]
		names(table)<-c('id',samplesnames1[ns])
		all<-merge(all,table,by="id",all=TRUE)
		
	}
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1
table1<-all
cutoffsamples1<-(length(samplesnames1)-1)*outl;

# print ((table1[1:10,]))
tablepvalue<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
tableratio<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
if (length(samplesnames1)>1){
	colpvalue<-1
# 	listloci<-vector("list",length=nrow(table))
	
	for (s1 in 1:(length(samplesnames1)-1)){
		for (s2 in (s1+1):length(samplesnames1)){
			ratio<-max1[s1]/max1[s2]
			tablepvalue[,colpvalue]<-mapply(normtest,((all[,s1+1])),((all[,s2+1])),ratio)
			tableratio[,colpvalue]<-(all[,s1+1])/(all[,s2+1])
			colpvalue<-colpvalue+1
		}
	}
	#do coherent conditions
	if (colpvalue>1){
		print(outl)
		print(cutoffsamples1)
		tempratio<-apply(tableratio,1,median)
		temppvalue<-apply(tablepvalue,1,median)
		numpval<-apply(tablepvalue,1,numpvalfunc)
		numexp<-apply(table1[2:(length(samplesnames1)+1)],1,numexpfunc)
		
		table1$tempr<-tempratio
		table1$tempp<-temppvalue
		table1$tempnp<-numpval
		table1$tempnexp<-numexp
		write.table(table1,"temp/group1",sep="\t",quote=F,row.names=F)
#print(table1[table1[,1]==1077,])
#		print(tablepvalue[table1[,1]==1077,])
		table1all<-table1[,c(1,ncol(table1))]
		names(table1all)<-c('id','t1')
		
 		table1<-table1[table1$tempnp<=cutoffsamples1 & table1$tempnexp<=outl,]
#		

#		print(table1[table1[,1]==1,])
#		print(table1[1:2,])
#
#print(table1[table1$id==280,])
		
		table1[,2]<-round(apply(table1[,2:(length(samplesnames1)+1)],1,median))
# 		print(table1[1:10,])
		table1<-table1[,1:2]
		
	}

}else{

	#temp<-apply(table[,2:(length(samplesnames1)+1)],1,mean)
	table1<-all
}
names(table1)<-c("id","freq1")
# table1[1:10,]
#  print ((table1))
# scale
#Load samples inside groups
# q()
ns<-0
listsamples2<-vector("list",length=length(samplesnames2))
max2<-1:length(samplesnames2)
#Load samples inside groups
for (s in samplesnames2){
	#s<-paste(sep="",s,projectval)
	ns<-ns+1
	print(s)
	infogroup<-append(infogroup,paste(sep="","G2:",s))

	query<-paste(sep="","select `idu`,`freq` from `",projectval,"`.`",s,"clusmap` where `idu`>0 group by `idu`;")
#  	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
	query<-paste(sep="","select SUM(freq) AS freq from `",projectval,"`.`",s,"clusraw` where `id`>0;")
#  	print(query)
	rs <- dbSendQuery(con,query) 
	max2[ns]<- unlist(as.vector(fetch(rs)))
# 	if (is.null(list$ref)==F){
# 		temp<-getcommon(temp)
# 	}
  		print(temp[1:5,])
		table<-applynorm(temp[,c(1,2)],scaleval,q1val,q2val,typetrval,typetdval)
# 		
		table<-merge(temp,table,by=1,all=FALSE)
		
#  	print(table[1:10,])
# 	print(listseqnames[1:2,])

		
	if (ns==1){
		all<-table[,c(1,3)]
		names(all)<-c('id',samplesnames2[ns])
# 		print(all[1:10,])
		#tableall<-temp
		
	}else{
		table<-table[,c(1,3)]
		names(table)<-c('id',samplesnames2[ns])
		all<-merge(all,table,by="id",all=TRUE)
		
	}
}

#calculate pvalue intra groups
all[is.na(all)]<-1
table2<-all
cutoffsamples2<-(length(samplesnames2)-1)*outl;

# print ((table2))
tablepvalue<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
tableratio<-matrix(ncol=sum(1:(ncol(all)-2)),nrow=nrow(all))
if (length(samplesnames2)>1){
	colpvalue<-1

	
	for (s1 in 1:(length(samplesnames2)-1)){
		for (s2 in (s1+1):length(samplesnames2)){
			ratio<-max2[s1]/max2[s2]
#			print ("statistic")
			tablepvalue[,colpvalue]<-mapply(normtest,((all[,s1+1])),((all[,s2+1])),ratio)
#print("norm")
			tableratio[,colpvalue]<-(all[,s1+1])/(all[,s2+1])
			colpvalue<-colpvalue+1
		}
	}
	
	#do coherent conditions
	if (colpvalue>1){
		
		tempratio<-apply(tableratio,1,median)
		temppvalue<-apply(tablepvalue,1,median)
		numpval<-apply(tablepvalue,1,numpvalfunc)
		numpexp<-apply(table2[2:(length(samplesnames2)+1)],1,numexpfunc)
		table2$tempr<-tempratio
		table2$tempp<-temppvalue
		table2$tempnp<-numpval
		table2$tempnexp<-numpexp
		
		table2all<-table2[,c(1,ncol(table2))]
		names(table2all)<-c('id','t2')
		write.table(table2,"temp/group2",sep="\t",quote=F,row.names=F)

#		print(table2[1:10,])
#table2<-table2[(table2$tempr>=cs & table2$tempr<=cb) & table2$tempnp>=cutoffsamples2 ,]
		table2<-table2[table2$tempnp<=cutoffsamples2 & table2$tempnexp<=outl,]

		table2[,2]<-round(apply(table2[,2:(ncol(table2)-2)],1,median))
#print(table2[1:10,])
		table2<-table2[,1:2]
		
	}
}else{

	#temp<-apply(table[,2:(length(samplesnames1)+1)],1,mean)
	table2<-all
}

names(table2)<-c("id","freq2")
###########remove not in sref
#  print (table1[1:10,])
#  print (table2[1:10,])
#pvalue of the two groups
table1[,3]<-table1[,2]
table2[,3]<-table2[,2]
table<-merge(table1,table2,by='id',all=TRUE)
table[is.na(table)]<-0

#print (table[1:10,])
checktable<-merge(table,table1all,by='id',all=TRUE)
checktable<-merge(checktable,table2all,by='id',all=TRUE)
checktable[is.na(checktable)]<-10-11	
#print (checktable[1:10,])
#write.table(checktable,"checktable",sep="\t",quote=F,row.names=F)

checktable2<-checktable[(checktable[,2]==0 & checktable$t1<0 & checktable[,4]>0) | (checktable[,4]==0 & checktable$t2<0 & checktable[,2]>0) | (checktable[,2]>0 & checktable[,4]>0),]
#print (checktable2[1:10,])
table<-checktable2[,1:5]

ratio<-mean(max1)/mean(max2)

if (typemetval=="Bayesian"){
	table$p<-mapply(statistic,table[,3],table[,5],ratio)
}else if (typemetval=="Binomial"){
	table$p<-mapply(binomial,table[,3],table[,5],mean(max1),mean(max2))
}else if (typemetval=="Ztest"){
	table$p<-mapply(ztest,table[,3],table[,5],sum(table[,3]),sum(table[,5]))
}else if (typemetval=="Fisher"){
	table$p<-mapply(fisher,table[,3],table[,5],sum(table[,3]),sum(table[,5]))
}else{
	table$p<-1
}
table$q<--1

if (corval!="0"){
	sort<-sort(table$p,index.return=T)
	ind<-1:nrow(table)
	order<-unlist(lapply(ind,function (x) ind[sort$ix==x]))
	table$q<-mapply(BHcorrection,table$p,order,nrow(table))
}else{
	table$q<-1
}
# # print (order)
# # print (sort)
# # print(table[1:10,])
# samplesnames2
# samplesnames1
table$ratio<-table[,3]/table[,5]
# table[is.infinite(table)==T]
table$ratio[is.infinite(table$ratio)==T]<-9999
table$ratio[is.na(table$ratio)]<-1
# q()
# des<-nrow(table[table$q<0.01 & (table$ratio>=1.5 | table$ratio<=0.5),])
# nrow(table)
# des
seqnames<-listseqnames
write.table(table,"table",sep="\t",quote=F,row.names=F)
#print(table)
table<-table[,c(1,2,4,6,7,8)]
namescols<-c("id",paste(collapse="-",samplesnames1),paste(collapse="-",samplesnames2),"pvalue","qvalue","ratio")
names(table)<-namescols
# write.table(file=paste(sep="",pathval,"table.txt"),table,col.names=T,row.names=F,quote=F,sep="\t")
# print(table[1:10,])


createtable<-paste(sep="","CREATE TABLE `difexp` (`id` int unsigned,`group 1` DOUBLE(8,2),`group 2` DOUBLE(8,2),`ratio` DOUBLE(8,2),pvalue DOUBLE(8,2),qvalue DOUBLE(8,2),INDEX (id));")
rs <- dbSendQuery(con,createtable) 

for (r in 1:nrow(table)){
 	
	
	loadtable<-paste(sep="","INSERT INTO `difexp` VALUES (",table[r,1],",",table[r,2],",",table[r,3],",",table$ratio[r],",",table$p[r],",",table$q[r],");")
#  	print (loadtable)
	rs <- dbSendQuery(con,loadtable) 
	
}

# 	
# 	write.table(table,file=name)
d<-getwd()
write.table(table,paste(sep="",d,"/","temp/df.done"))
dbDisconnect(con) 

}
d<-getwd()
name<-paste(sep="",d,"/","temp/args.txt")
d<-read.table(name)
# d[,1]
opt<-as.vector(d[,1])
opt
g1<-unlist(strsplit(opt[6],":"))
g1
g2<-unlist(strsplit(opt[7],":"))
g2
source("Rscripts/R/db.R")
# Ztest,BH,1000000,0,100,test2:,test2:,testpro,scale,localhost,lpantano,sqllorena
 difexpseq(opt[1],opt[2],opt[3],0,100,g1,g2,opt[8],"",opt[9],hostname,username,pssw,port,opt[10])
# difexpseq("Ztest","BH","1000000",0,100,"test2","test2","testpro","","scale","localhost","lpantano","sqllorena")
q()
