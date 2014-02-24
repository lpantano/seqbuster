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
Do DE analysis with DESeq
}

