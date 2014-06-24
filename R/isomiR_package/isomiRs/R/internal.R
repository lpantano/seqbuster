#' put header to input files
#'
#' @param table table from IsomirDataSeq
put.header<-function(table)
{
    names(table)[c(1,3,4,7,8,9,10,13,14)]<-c("seq","freq","mir","mism","add","t3","t5","DB","ambiguity")
    table<-table[,c(1,3,4,7,8,9,10,13,14)]
    table[,2]<-as.numeric(table[,2])
    return(table)
}

#' filter by relative abundance to reference
#' @import data.table
#' @param table object miraligner table
#' @param limit remove sequences lower than this number
filter.by.cov<-function(table,limit=10)
{
    freq=NULL
    tab.fil<-table[table$DB=="miRNA",]
    tab.fil.out<-data.table(tab.fil)
    tab.fil.out<-as.data.frame(tab.fil.out[,list(total=sum(freq)),
        by=c('mir')])
    tab.fil<-merge(tab.fil[,c(3,1:2,4:ncol(tab.fil))],tab.fil.out,
        by=1)
    tab.fil$score<-tab.fil$freq/tab.fil$total*100
    tab.fil<-tab.fil[tab.fil$score>=limit,]
    return (tab.fil)
}

#' Filter tablo reference
#'
#' @param table object miraligner table
#' @param cov remove sequences that have relative abundance lower than this number
filter.table<-function(table,cov=10)
{
    table<-put.header(table)
    table<-filter.by.cov(table,cov)
    return(table)
}


#' plot general information
#'
#' @param table object miraligner table
#' @param colid column to be analize
isomir.general.type<-function(table,colid)
{
    temp<-table
    temp$idfeat<-paste(table[,colid],table$mir)
    temp<-temp[order(temp$idfeat),]
    temp<-temp[!duplicated(temp$idfeat),]
    temp<-as.data.frame(summary(temp$mir))
    feat.dist<-cut(as.numeric(temp[,1]),breaks=c(-1,0.5,1.5,2.5,Inf),labels=c("0","1","2",">3"))
    return (as.data.frame(summary(feat.dist)))  
}

#' do counts table considering what isomiRs take into account
#'
#' @param x object isomiDataSeq
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param mism differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
do.mir.table<-function(x,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F)
{
    table.merge<-data.frame()
    des<-x@design
    for (sample in row.names(des)){
        print (sample)
        d<-x@expList[[sample]]
        d<-collapse.mirs(d,ref=ref,iso5=iso5,iso3=iso3,add=add,
            mism=mism,seed=seed)
        names(d)[ncol(d)]<-sample 
        if( nrow(table.merge)==0){
            table.merge<-d
        }else{
            table.merge<-merge(table.merge,d,by=1,all=T)
        }        
    }
    row.names(table.merge)<-table.merge[,1]
    table.merge<-as.matrix(table.merge[,2:ncol(table.merge)])
    table.merge[is.na(table.merge)]<-0
    x@counts<-as.matrix(table.merge)
    return(x)
}

#' Collapse isomiRs in miRNAs 
#' @import data.table
#' @param table object miraligner table
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param mism differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
collapse.mirs<-function(table,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F)
{
    label<-table$mir
    freq=NULL
    if (ref==T){
        ref.val<-do.call(paste,table[,4:7])
        ref.val[grep("[ATGC]",ref.val,invert=T)]<-"ref"
        ref.val[grep("[ATGC]",ref.val)]<-"iso"
        label<-paste(label,ref.val,sep=".")
    }
    if (iso5==T){
        label<-paste(label,table[,7],sep=".")
    }
    if (seed==T){
        seed.val<-as.character(table[,4])
        seed.val[grep("^[2-8][ATGC]",seed.val,invert=T)]<-"no"
        label<-paste(label,seed.val,sep=".")
    }
    if (iso3==T){
        label<-paste(label,table[,6],sep=".")
    }
    if (add==T){
        label<-paste(label,table[,5],sep=".")
    }
    if (mism==T){
        label<-paste(label,table[,4],sep=".")
    }
    
    table$id<-label
    table.out<-data.table(table)
    table.out<-as.data.frame(table.out[,list(total=sum(freq)),by=c('id')])
    table.out[is.na(table.out)]<-0
    return(table.out)
}

#' Do summary of different isomiRs events
#'
#' @param table object miraligner table
#' @param colid column to be considered
isomir.position<-function(table,colid)
{
    temp<-table
    temp[,colid]<-as.character(temp[,colid])
    pos<-as.data.frame(t(as.data.frame((strsplit(temp[,colid],"-",fixed=2)))))
    row.names(pos)<-1:nrow(pos)
    pos$mir<-temp$mir
    pos$freq<-temp$freq
    pos<-pos[pos[,1]!=0,]
    pos$size<-apply(pos,1,function(x){
        p<-length(unlist(strsplit(x[2],"")))
        if(x[1]=="d"){
            p<-p*-1
        }
        return(p)
    })
    pos$idfeat<-paste(pos$size,pos$mir)
    pos<-pos[order(pos$idfeat,abs(pos$size)),]
    return (pos[,c(3,4,5)])  
}

#' Do summary of nt substitution events
#'
#' @param table object miraligner table
#' @param colid column to be considered
subs.position<-function(table,colid)
{
    temp<-table
    temp[,colid]<-as.character(temp[,colid])
    nt<-sub("[0-9]+","",temp[,colid])  
    pos<-sub("[ATGC]{2}","",temp[,colid])
    pos<-data.frame(nt=as.character(nt),size=pos,mir=temp$mir,freq=temp$freq)
    pos$nt<-as.character(pos$nt)
    pos<-pos[pos[,1]!="",]
    nt.2<-as.data.frame(t(as.data.frame((strsplit(pos$nt,"",fixed=2)))))
    names(nt.2)<-c("current","reference")
    names(nt.2$current)<-""
    names(nt.2$reference)<-""
    pos<-cbind(pos,nt.2)
    pos$size<-factor(pos$size,levels=1:25)
    return (pos[,c(3,4,2,5,6)]) 
}
