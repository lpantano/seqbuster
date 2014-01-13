

parse<-function (txt)
{

#read file
param<-read.table(txt,sep="\t")
name<-vector()
list<-vector("list",length=nrow(param))
for (i in 1:nrow(param)){
	opt<-unlist(strsplit(as.character(param[i,1]),":"))	
	name<-append(name,opt[1])
	if (length(grep("group",opt[1]))>0){
		list[[i]]<-unlist(strsplit(opt[2],","))
	}else{
		list[[i]]<-opt[2]
	}	

}

names(list)<-name
return(list)

}

