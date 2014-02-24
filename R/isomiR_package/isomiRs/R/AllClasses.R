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
#' create class
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
    listObj[[row.names(design)[idx]]]<-d[,c(1:3,6:9)]
    listObjVar[[row.names(design)[idx]]]<-out
  }
  IsoObj@design<-design
  IsoObj@expList<-listObj
  IsoObj@varList<-listObjVar
  IsoObj<-do.mir.table(IsoObj)
  #class(listObj)<-"Isomirs"
  return(IsoObj)
}
