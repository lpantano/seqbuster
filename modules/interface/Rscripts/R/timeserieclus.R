
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
	print(freq)
# 	print(nrow(data))
	data<-data[log(data$freq,base=2)>=freq,]
# 	print(nrow(data))
	return (data)
}

applynorm<-function(tempn,scaleval,q2val,q1val,typetrval,typetdval){
#  	print(c(scaleval,q1val,q2val,typetrval,typetdval))
	
	if (as.numeric(q1val)<100 | as.numeric(q2val)>0 ){
		tempn<-discartfreq(tempn,as.numeric(q1val),as.numeric(q2val))
# 		q<-quantile(tempn[,2],probs=c(10/100,90/100))
# 		print(tempn[1:10,])
		#print("quantile")
# 		print(tempn[1:10,])
	}
	
	if(typetdval!=""){
		#print("center")
#  		print(tempn[1:10,])
		tempn<-centernorm(tempn,typetdval)
#  		print(tempn[1:10,])
	}
	if (scaleval!="" ){
# 		print(scaleval)
#  		print(tempn[1:10,])
		tempn$freq<-scalefreq(tempn$freq,as.numeric(scaleval))
# 		print(tempn[1:10,])
	}
	if (typetrval!=""){
		#print("trans")
		tempn$freq<-logfreq(tempn$freq,(typetrval))
# 		print(tempn[1:10,])
	}
	
	

	
	return(tempn)
}

LmListInt1<-function (y, par, m) 
   { lm(y ~ I(par[[1]] * m[1]) + I(par[[2]] * m[2]) + I(par[[3]] * m[3])) }
library(RMySQL)

# Ztest,BH,1000000,0,100,test2:,test2:,testpro,scale,localhost,lpantano,sqllorena
difexpseq<-function(scaleval,q1val,q2val,samplesnames1,projectval,typetrval,typetval,localhost,user,pssw,val){
# print (c(typemetval,corval,scaleval,q1val,q2val,samplesnames1,samplesnames2,projectval,typetrval,typetdval,localhost,user))
MySQL(max.con = 16, fetch.default.rec = 5000000, force.reload = TRUE)

m <- dbDriver("MySQL")

con <- dbConnect(m,host=localhost,user=user,db=projectval,password=pssw)


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
	#print(s)
	infogroup<-append(infogroup,paste(sep="","G1:",s))

	query<-paste(sep="","select `idu`,`freq` from `",projectval,"`.`",s,"clusmap` where `freq`>0 and `idu`>0 group by `chr`;")
#  	print(query)
	rs <- dbSendQuery(con,query) 
	temp <- as.data.frame(fetch(rs))
# 	if (is.null(list$ref)==F){
# 		temp<-getcommon(temp)
# 	}
#   		print(temp[1:10,c(1,2)])
		table<-applynorm(temp,scaleval,q1val,q2val,typetrval,"")
		if (nrow(table)==0){
			table<-data.frame(idu=0,freq=0)
			names(table)<-c('idu',s)
		}

# 		print(table[1:10,c(1,3)])
if (ns==1){
		all<-table
		names(all)<-c('idu',s)
# 		print(all[1:10,])
		#tableall<-temp
		
	}else{
		
		names(table)<-c('idu',s)
		all<-merge(all,table,by="idu",all=TRUE)
		
	}
}
# q()
#calculate pvalue intra groups
all[is.na(all)]<-1

parI = list(val, val*val, val*val*val)

Mat<-(all[,2:ncol(all)])

Genes = nrow(Mat); PI = length(parI); LI = 2^ (PI); 
Models = matrix(0, 8, 3)
{ j = 1     
{
for (i1 in 0:1) {
   for (i2 in 0:1) {
      for (i3 in 0:1) {
         Models[j,] = c(i1,i2,i3)
               j = j + 1 } } } } }
PCUTOFF<-0.05
Result = t(apply(Mat, 1, function(y) {
         lr = lapply(1:8, function(j) {LmListInt1(y, parI, Models[j,])})
         mod = sapply(2:8, function(k) {anova(lr[[1]], lr[[k]])$Pr[2]})
         mod[is.na(mod)] = 1
         c(ifelse(min(mod) <= PCUTOFF, which.min(mod)+1, 1), min(mod))
      }))
      rownames(Result) = rownames(Mat)
      colnames(Result) = c("Model", "P-value")

table<-cbind(all,Result)

# write.table(midfr,file="/projects/srna_sps/affy/time.serie.data.age",row.names=F,quote=F,sep="\t")

# table$type<-mapply(getid,table[,1])
# print(table[1:10,])
# table<-table[,c(1,2,4,6,7,8)]
# namescols<-c("id",paste(collapse="-",samplesnames1),"model","pvalue")
# names(table)<-namescols
# write.table(file=paste(sep="",pathval,"table.txt"),table,col.names=T,row.names=F,quote=F,sep="\t")
# print(table[1:10,])
#print(table[1,])
 query<-paste(sep="","DROP TABLE IF EXISTS `",projectval,"`.`timeserie`;")
 rs <- dbSendQuery(con,query) 

createtable<-paste(sep="","CREATE TABLE `",projectval,"`.`timeserie` (`id` int unsigned, ")
for (c in 2:(ncol(table)-2)){
	createtable<-paste(sep="",createtable,"`",names(table)[c],"` DOUBLE(5,2), ")
}
createtable<-paste(sep="",createtable," `model` int,pvalue DOUBLE(8,2),INDEX (id));")
#print (createtable)
 rs <- dbSendQuery(con,createtable) 

for (r in 1:nrow(table)){
 	
	
	loadtable<-paste(sep="","INSERT INTO `",projectval,"`.`timeserie` VALUES (")
	for (c in 1:(ncol(table)-1)){
			loadtable<-paste(sep="",loadtable,table[r,c],", ")
	}
	loadtable<-paste(sep="",loadtable,table[r,ncol(table)]," );")

   	#print (loadtable)
 	rs <- dbSendQuery(con,loadtable) 
	
}

# 	
d<-getwd()
write.table(table,paste(sep="",d,"/","temp/ts.done"))
dbDisconnect(con) 

}
d<-getwd()
name<-paste(sep="",d,"/","temp/args.txt")
d<-read.table(name)
# d[,1]
opt<-as.vector(d[,1])
# opt
g1<-unlist(strsplit(opt[1],":"))
#g1<-unlist(strsplit("hsa2dspsmir:hsa204dspsmir:hsa5105dspsmir:hsa9277dspsmir",":"))

#  difexpseq(opt[1],opt[2],opt[3],0,100,g1,g2,opt[8],"",opt[9],opt[10],opt[11],opt[12])
#  scaleval,q1val,q2val,samplesnames1,projectval,typetrval,typetdval,localhost,user,pssw
# vec<-c(2,204,5101,9277)
vec<-as.numeric(unlist(strsplit(opt[6],",")))

difexpseq(1000000,0,100,g1,opt[2],"log2","scale",opt[3],opt[4],opt[5],vec)
 

