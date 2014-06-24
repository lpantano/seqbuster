#' Do DE analysis with DESeq
#' 
#' This function does differential expression analysis with DESeq2.
#' It will return a DESeq2 object.
#' You can merge all isomiRs into miRNA by calling the function only 
#' with the frist two parameters. You can get a table with isomiRs and 
#' the reference miRBase sequence by calling the function with ref=T. 
#' You can get a table with 5' trimming isomiRS, miRBase reference and 
#' the rest by calling with ref=T,iso5=T. 
#' If you set up all parameters to TRUE, you will get a table for 
#' each differnt sequence mapping to a miRNA
#' 
#' @param x object isomiDataSeq
#' @param formula used for DE analysis
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param mism differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
#' @return DESeq object
#' @examples
#' library(DESeq2)
#' data(isomiRex)
#' dds<-deIso(isomiRex,formula=~condition)
#' @export
#' @import DESeq2
deIso<-function(x,formula,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F)
{
    if (ref | iso5 | iso3 | add | mism | seed){
        x<-do.mir.table(x,ref,iso5,iso3,add,mism,seed)
    }
    countData<-x@counts
    dds<-DESeqDataSetFromMatrix(countData = countData,
        colData = x@design,
        design = formula)
    dds <- DESeq(dds,quiet=T)
    dds
}

#' Plot the amount of isomiRs in different samples divided by group
#'
#' @param x object isomiDataSeq
#' @param top number of isomiRs used
#' @export
#' @import ggplot2
#' @import gplots
#' @import RColorBrewer
plotTop<-function(x,top=20)
{
    dds<-x[["dds"]]
    rld <- rlogTransformation(dds)
    #vsd <- varianceStabilizingTransformation(dds)
    select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:top]
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(counts(dds,normalized=TRUE), col = hmcol,
              scale="none",
              dendrogram="none", trace="none")
}

#' Plot the amount of isomiRs in different samples divided by group
#'
#' @param x object isomirDataSeq
#' @param type string (t5,t3,add,sub) to indicate what isomiR 
#' change to use for the plot
#' 
#' @return ggplot2 figure showing the selected isomiR changes among samples 
#' @examples
#' library(isomiRs)
#' data(isomiRex)
#' isomiRex<-plotIso(isomiRex,"t5")
#' @export
#' @import ggplot2
#' @examples
#' data(isomiRex)
#' plotIso(isomiRex)
plotIso<-function(x,type="t5")
{
    freq=size=group=abundance=NULL
    codevn<-c(2,3,4,5)
    names(codevn)<-c("t5","t3","sub","add")
    ratiov<-c(1/6,1/6,1/23,1/3)
    names(ratiov)<-names(codevn)
    coden<-codevn[type]
    ratio<-ratiov[type]
    des<-x@design
    table<-data.frame()
    for (sample in row.names(des)){
        print(sample)
        temp<-data.table(x@varList[[sample]][[coden]])
        uniq.dat<-as.data.frame(table(x@varList[[sample]][[coden]]$size))
        temp<-as.data.frame(temp[,list(freq=sum(freq)),by="size"])
        total<-sum(temp$freq)
        temp<-merge(temp,uniq.dat,by=1)
        Total<-sum(temp$Freq)
        temp$abundance<-temp$freq/total
        temp$unique<-temp$Freq/Total
        table<-rbind(table,data.frame(size=temp$size,abundance=temp$abundance,
            unique=temp$unique,sample=rep(sample,nrow(temp)),
            group=rep(des[sample,"condition"],nrow(temp))))
    }
    x@sumList[[type]]<-table
    p <- ggplot(table)+
        geom_jitter(aes(x=factor(size),y=unique,colour=factor(group),
            size=abundance))+
        scale_colour_brewer("Groups",palette="Set1")+
        theme_bw(base_size = 14, base_family = "") +
        theme(strip.background=element_rect(fill="slategray3"))+
        labs(list(title=paste(type,"distribution"),y="# of unique sequences",
            x="position respect to the reference"))
    print(p)  
    x
}

#' create count tables from isomirs
#' @param x IsomirDataSeq class
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param mism differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
#' @return count table
#' 
#' @export
makeCounts<-function(x,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F)
{
    x<-do.mir.table(x,ref,iso5,iso3,add,mism,seed)
    x
}



#' normalize count data
#'
#' @param x IsomirDataSeq object
#' @param formula formula that will be used for DE
#' @export
normIso<-function(x,formula=~condition)
{
    countData<-x@counts
    dds<-DESeqDataSetFromMatrix(countData = countData,
                                colData = x@design,
                                design = formula)
    rld<-rlogTransformation(dds,blind=FALSE)
    x@normcounts <- assay(rld)
    x
}

