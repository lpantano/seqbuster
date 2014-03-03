#' IsomirDataSeq object and constructors
#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq<-setClass("IsomirDataSeq",
         slots=c(counts="matrix",
                 design="data.frame",
                 varList="list",
                 expList="list",
                 sumList="list"
           ))
#' load data
#'
#' @param files all samples
#' @param cov remove sequences that have relative abundance lower than this number
#' @param design data frame containing groups for each sample
#' @param header files contain headers
#' @param skip skip first line when reading files
#' 
#' @export
loadIso2<-function(files,design,cov=1,header=F,skip=1){
  IsoObj<-IsomirDataSeq()
  listObj<-vector("list")
  listObjVar<-vector("list")
  idx<-0
  for (f in files){
    idx<-idx+1
    print(idx)
    d<-read.table(f,header=header,skip=skip)
    
    d<-filter.table(d,cov)
    #Run function to describe isomirs
    out<-list(summary=0,t5sum=isomir.position(d,6),t3sum=isomir.position(d,7),subsum=subs.position(d,4),addsum=isomir.position(d,5))
    #class(out)<-"Isomirs"
    listObj[[row.names(design)[idx]]]<-d
    listObjVar[[row.names(design)[idx]]]<-out
  }
  IsoObj@design<-design
  IsoObj@expList<-listObj
  IsoObj@varList<-listObjVar
  IsoObj<-do.mir.table(IsoObj)
  #class(listObj)<-"Isomirs"
  return(IsoObj)
}

#' create count tables from isomirs
#'
#' @param IsomirDataSeq class
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param mism differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
#' @return count table
#' 
#' @export
makeCounts<-function(x,ref=F,iso5=F,iso3=F,add=F,mism=F,seed=F){
  x<-do.mir.table(x,ref,iso5,iso3,add,mism,seed)
  return(x@counts)
}



