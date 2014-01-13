
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

bootstrapin<-function(n1,n2,m1,m2){



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
# print(c(nmv,mv))
# print(nrow(data))
# print((data))
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
	print (sum(freqs))
	newfreqs<-round(freqs/sum(freqs)*sc)
	return (newfreqs)
}


logfreq<-function(freqs,l){
	if (l=="n"){
		l=exp(1)
	}
	if (l=="log2"){
		l=2
	}
	if (l=="log10"){
		l=10
	}
	
	l<-as.numeric(l)
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
	print(freq)
# 	print(nrow(data))
	data<-data[log(data$freq,base=2)>=freq,]
# 	print(nrow(data))
	return (data)
}

applynorm<-function(tempn,list){
# 	print(tempn)
	
	if (as.numeric(list$qmax)<100 | as.numeric(list$qmin)>0 ){
		tempn<-discartfreq(tempn,as.numeric(list$qmax),as.numeric(list$qmin))
# 		q<-quantile(tempn[,2],probs=c(10/100,90/100))
	
	}
	
	if(is.null(list$center)==F){
# 		print(tempn[1:10,])
		tempn<-centernorm(tempn,list$center)
# 		print(tempn[1:10,])
	}
	if (is.null(list$scale)==F){
		
		tempn$freq<-scalefreq(tempn$freq,as.numeric(list$scale))
		
	}
	if (is.null(list$trans)==F){
# 		print(tempn[1:10,])
		tempn$freq<-logfreq(tempn$freq,(list$trans))
# 		print(tempn[1:10,])
	}
	
	
# 	if (list$norm=="hyper"){
# 		tempn<-normbyhyper(tempn,as.numeric(list$scale),as.numeric(list$qmax),as.numeric(list$qmin))
# # 		print(tempn[1:10,])
# 	}
	
	
	
	
	
	return(tempn)
}


