library(data.table)


do.mir.table<-function(files,names,cov=10,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F){
  table.merge<-data.frame()
  for (f in files){
    
     d<-read.table(f,header=T)
     d<-put.header(d)
     d<-filter.by.cov(d,cov)
     d<-collapse.mirs(d,ref=ref,iso5=iso5,iso3=iso3,add=add,mism=mism,seed=seed)
     names(d)[ncol(d)]<-names[files==f]
      
    if( nrow(table.merge)==0){
      table.merge<-d
    }else{
      table.merge<-merge(table.merge,d,by=1,all=T)
    }
    
  }
  return(table.merge)
}


put.header<-function(table){
  
  if (names(table)[3]!="mir"){
    names(table)[c(2,3,6,7,8,9,12,13)]<-c("freq","mir","mism","add","t3","t5","DB","ambiguity")
  }
  return(table)
}

filter.by.cov<-function(table,limit=10){
  
  table.out<-data.table(table)
  table.out<-table.out[table.out$DB=="miRNA",]
  table.out<-as.data.frame(table.out[,list(total=sum(freq)),by=c('mir')])
  table<-merge(table[,c(3,1:2,4:ncol(table))],table.out,by=1)
  table$score<-table$freq/table$total*100
  table<-table[table$score>=limit,]
  return (table)
  
}


collapse.mirs<-function(table,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F){
  label<-table$mir
  if (ref==T){
    ref.val<-do.call(paste,table[,6:9])
    ref.val[grep("[ATGC]",ref.val,invert=T)]<-"ref"
    ref.val[grep("[ATGC]",ref.val)]<-"iso"
    label<-paste(label,ref.val,sep=".")
  }
  if (iso5==T){
    label<-paste(label,table[,8],sep=".")
  }
  if (seed==T){
    seed.val<-as.character(table[,6])
    seed.val[grep("^[2-8][ATGC]",seed.val,invert=T)]<-"no"
    label<-paste(label,seed.val,sep=".")
  }
  if (iso3==T){
    label<-paste(label,table[,9],sep=".")
  }
  if (add==T){
    label<-paste(label,table[,7],sep=".")
  }
  if (mism==T){
    label<-paste(label,table[,6],sep=".")
  }
  
  table$id<-label
  table.out<-data.table(table)
  
  table.out<-as.data.frame(table.out[,list(total=sum(freq)),by=c('id')])
  table.out[is.na(table.out)]<-0
  return(table.out)
}
