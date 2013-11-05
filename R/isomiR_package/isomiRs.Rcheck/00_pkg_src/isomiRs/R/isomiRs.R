library(data.table)
library(ggplot2)


do.mir.table<-function(files,names,cov=10,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F){
  table.merge<-data.frame()
  for (f in files){
    
    d<-read.table(f,header=T)
    d[,2]<-as.numeric(d[,2])
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

isomir.general.type<-function(table,colid){
  temp<-table
  temp$idfeat<-paste(table[,colid],table$mir)
  temp<-temp[order(temp$idfeat),]
  temp<-table[!duplicated(temp$idfeat),]
  temp<-as.data.frame(summary(temp$mir))
  feat.dist<-cut(as.numeric(temp[,1]),breaks=c(-1,0.5,1.5,2.5,Inf),labels=c("0","1","2",">3"))
  return (as.data.frame(summary(feat.dist)))  
}


make.figure.general<-function(table){
  
  mism<-isomir.general.type(table,6)
  add<-isomir.general.type(table,7)
  t3<-isomir.general.type(table,8)
  t5<-isomir.general.type(table,9)
  dat<-data.frame(rbind(t5,t3,add,mism))
  names(dat)<-"count"
  dat$type<-c(rep("5 trimming",nrow(t5)),rep("3 trimming",nrow(t3)),rep("nt addition",nrow(add)),rep("nt subst",nrow(mism)))
  dat$iso<-c(row.names(t5),row.names(t3),row.names(add),row.names(mism))
  dat$iso<-factor(dat$iso,levels=c(0,1,2,">3"))
  ggplot(dat,aes(x = factor(type), y = count, fill=iso)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette="Set1") +
    theme_bw(base_size = 12, base_family = "") +
    labs(list(x="",y="# of miRNAs",fill="# of isomirs"))
}

filter.table<-function(table,cov=10){
  
  table[,2]<-as.numeric(table[,2])
  table<-put.header(table)
  table<-filter.by.cov(table,cov)
  return(table)
  
}