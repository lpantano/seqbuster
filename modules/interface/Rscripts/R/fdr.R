
library("RMySQL")
source("Rscripts/R/parse.R")
source("Rscripts/R/stand.R")

source("Rscripts/R/db.R")

d<-getwd()
name<-paste(sep="",d,"/","temp/usRNA")
d<-read.table(name)
opt<-as.vector(d[,1])
project<-opt[1]
sample<-opt[2]

MySQL(max.con = 16, fetch.default.rec = 5000000, force.reload = TRUE)

m <- dbDriver("MySQL")

if (port==0){
	con <- dbConnect(m,host=hostname,user=username,db=project,password=pssw)
}else{
	con <- dbConnect(m,host=hostname,user=username,db=project,password=pssw,port=port)
}

query<-paste(sep="","select `id`,`freq` from `",project,"`.`",sample,"clusraw` where `freq` > 0;")
rs <- dbSendQuery(con,query) 
total <- as.data.frame(fetch(rs))

deeph<-sum(total[,2])

query<-paste(sep="","select `id`,`freq` from `",project,"`.`",sample,"clusmap` where `id`>0 group by `id`;")
rs <- dbSendQuery(con,query) 
temp <- as.data.frame(fetch(rs))
fdr<-temp

for (i in 1:nrow(temp)){
 counts<-temp[i,2]
 prob<-counts/2/deeph
 currentd<-sum(total[total$freq > counts,2])
 coverage<-round(deeph-currentd,digits=0)
 if (coverage>counts){
  pval<-1-sum(dbinom(1:counts,coverage,prob))
 }else{
   pval<-1	
 }
 pval<-formatC(pval, format="e",digits=2)
 fdr[i,2]<-pval
}

write.table(fdr,"temp/fdr",row.names=F,col.names=F,quote=F,sep="\t")

dbDisconnect(con) 


