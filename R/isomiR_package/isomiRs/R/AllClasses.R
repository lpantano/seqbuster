IsomirDataSeq<-setClass("IsomirDataSeq",
         slots=c(counts="matrix",
                 design="data.frame",
                 varList="list",
                 expList="list",
                 sumList="list"
           ))

loadIso2<-function(files,design,cov=10,header=F,skip=1){
  IsoObj<-IsomirDataSeq()
  listObj<-vector("list")
  listObjVar<-vector("list")
  idx<-0
  for (f in files){
    idx<-idx+1
    print(idx)
    d<-read.table(f,header=header,skip=skip)
    d[,2]<-as.numeric(d[,2])
    d<-put.header(d)
    d<-filter.by.cov(d,cov)
    #Run function to describe isomirs
    out<-list(summary=0,t5sum=isomir.position(d,8),t3sum=isomir.position(d,9),subsum=subs.position(d,6),addsum=isomir.position(d,7))
    #class(out)<-"Isomirs"
    listObj[[row.names(design)[idx]]]<-d[,c(1:3,6:9)]
    listObjVar[[row.names(design)[idx]]]<-out
  }
  IsoObj@design<-design
  IsoObj@expList<-listObj
  IsoObj@varList<-listObjVar
  #class(listObj)<-"Isomirs"
  return(IsoObj)
}
