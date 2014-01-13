
source("Rscripts/R/parse.R")
list<-parse("Rscripts/script_param")
write.table(file="Rscripts/testoutput",list)
