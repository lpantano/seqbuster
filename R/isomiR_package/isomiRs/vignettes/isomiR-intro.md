<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{An Introduction to the IsomiR package}
-->

Example of isomiRs function 

========================================================

We are going to use a small RNAseq data coming from human frontal cortex samples to give some basic examples of isomiRs analyses. There are two groups: b) individuals with less than a year, o) individuals in the eldery.



```r
library(isomiRs)
```

```
## Error: there is no package called 'isomiRs'
```

```r
data(isomiRex)
```

```
## Warning: data set 'isomiRex' not found
```

You can plot isomiRs tendencies with `plotIso`


```r
#check plot iso
obj<-plotIso(isomiRex,type="t5")
```

```
## Error: could not find function "plotIso"
```

Do differential expression analysis using DESeq2


```r
dds<-deIso(obj,formula=~condition,ref=TRUE,iso5=T)
```

```
## Error: could not find function "deIso"
```

```r
#plotMA
library(DESeq2)
```

```
## Loading required package: GenomicRanges
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: IRanges
## Loading required package: XVector
## Loading required package: Rcpp
## Loading required package: RcppArmadillo
```

```r
plotMA(dds)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'plotMA': Error: object 'dds' not found
```

Let's do a PLS model having in the count matrix the reference miRNA, and all isomiRs for each miRNA.


```r
obj<-makeCounts(obj,ref=T)
obj<-normIso(obj)
pls.obj<-isoPLSDA(obj,"condition")
isoPLSDAplot(pls.obj$component,obj@design[,"condition"])
```

You can do the analysis just with features that are important for the model


```r
pls.obj<-isoPLSDA(obj,"condition",refinment=TRUE)
```
