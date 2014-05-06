\name{deIso}
\alias{deIso}
\title{Do DE analysis with DESeq}
\usage{
deIso(x, formula, ref = F, iso5 = F, iso3 = F, add = F, mism = F,
  seed = F)
}
\arguments{
  \item{x}{object isomiDataSeq}

  \item{formula}{used for DE analysis}

  \item{ref}{differenciate reference miRNA from rest}

  \item{iso5}{differenciate trimming at 5 miRNA from rest}

  \item{iso3}{differenciate trimming at 3 miRNA from rest}

  \item{add}{differenciate additions miRNA from rest}

  \item{mism}{differenciate nt substitution miRNA from
  rest}

  \item{seed}{differenciate changes in 2-7 nt from rest}
}
\value{
DESeq object
}
\description{
This function does differential expression analysis with
DESeq2. It will return a DESeq2 object. You can merge all
isomiRs into miRNA by calling the function only with the
frist two parameters. You can get a table with isomiRs and
the reference miRBase sequence by calling the function with
ref=T You can get a table with 5' trimming isomiRS, miRBase
reference and the rest by calling with ref=T,iso5=T If you
set up all parameters to TRUE, you will get a table for
each differnt sequence mapping to a miRNA
}
\examples{
library(DESeq2)
data(isomiRex)
dds<-deIso(isomiRex,formula=~condition)
}

