

summary<-function(x){
  UseMethod("summary")
}


plotTop<-function(x,top=20){
  UseMethod("plotTop")
}
lmsIso<-function(x,formula,merge="all"){
  UseMethod("summary")
}

coexpIso<-function(x,merge="all"){
  UseMethod("summary")
}


summary.Isomirs<-function(x){
  #whatever to do with my object (generic information)
}


deIso<-function(x,iso5=F,iso3=F,add=F,mism=F,seed=F,formula= ~condition){
  print("doing")
  countData<-do.mir.table(x,iso5=iso5,iso3=iso3,add=add,mism=mism,seed=seed)
  x[["countData"]]<-countData
  dds<-DESeqDataSetFromMatrix(countData = countData,
                         colData = x[["design"]],
                         design = formula)
  dds <- DESeq(dds,quiet=T)
  x[["dds"]]<-dds
  return(x)
}


plotTop.Isomirs<-function(x,top=20){
  dds<-x[["dds"]]
  rld <- rlogTransformation(dds)
  #vsd <- varianceStabilizingTransformation(dds)
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:top]
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(counts(dds,normalized=TRUE), col = hmcol,
  scale="none",
  dendrogram="none", trace="none")
  
  return(1)
}

plotIso<-function(x,type="t5"){
  #whatever to do with my object (generic information)
  codev<-c(4,5,6,7)
  names(codev)<-c("t5","t3","sub","add")
  ratiov<-c(1/6,1/6,1/23,1/3)
  names(ratiov)<-names(codev)
  code<-codev[type]
  ratio<-ratiov[type]
  des<-x[["design"]]
  table<-data.frame()
  #print(paste(code,ratio))
  for (sample in row.names(des)){
    print(sample)
    temp<-data.table(x[[sample]][[code]])
    uniq.dat<-as.data.frame(table(x[[sample]][[code]]$size))
    temp<-as.data.frame(temp[,list(freq=sum(freq)),by="size"])
    total<-sum(temp$freq)
    temp<-merge(temp,uniq.dat,by=1)
    Total<-sum(temp$Freq)
    #temp$freqn<-log2(temp$freq/(ratio*total))
    temp$abundance<-temp$freq/total
    temp$unique<-temp$Freq/Total
    table<-rbind(table,data.frame(size=temp$size,abundance=temp$abundance,
                                  unique=temp$unique,sample=rep(sample,nrow(temp)),
                                  group=rep(des[sample,"condition"],nrow(temp))))
  }
  x[[type]]<-table
  #print(table)
  p <- ggplot(table)+
    geom_jitter(aes(x=factor(size),y=unique,colour=factor(group),
                    size=abundance))+
    #geom_jitter(position=position_jitter(width=0.1),
    #            aes(factor(size),score,colour=sample)) +
    #scale_fill_brewer(palette="Set1",guide="none")+
    scale_colour_brewer(palette="Set1",guide="none")+
    theme_bw(base_size = 14, base_family = "") +
    theme(strip.background=element_rect(fill="slategray3"))+
    labs(list(title=paste(type,"distribution"),y="$ of unique sequences",
              x="position respect to the reference"))
    print(p)  
    return(x)
}


loadIso<-function(files,design,cov=10,header=F,skip=1){
  listObj<-vector("list")
  idx<-0
  for (f in files){
    idx<-idx+1
    print(idx)
    d<-read.table(f,header=header,skip=skip)
    d[,2]<-as.numeric(d[,2])
    d<-put.header(d)
    d<-filter.by.cov(d,cov)
    #Run function to describe isomirs
    out<-list(counts=d,design=design[idx,],summary=0,te5sum=isomir.position(d,8),t3sum=isomir.position(d,9),subsum=subs.position(d,6),addsum=isomir.position(d,7))
    #class(out)<-"Isomirs"
    listObj[[row.names(design)[idx]]]<-out
  }
  listObj[["design"]]<-design
  class(listObj)<-"Isomirs"
  return(listObj)
}
  
