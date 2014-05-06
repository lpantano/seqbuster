\name{makeCounts}
\alias{makeCounts}
\title{create count tables from isomirs}
\usage{
makeCounts(x, ref = F, iso5 = F, iso3 = F, add = F, mism = F,
  seed = F)
}
\arguments{
  \item{x}{IsomirDataSeq class}

  \item{ref}{differenciate reference miRNA from rest}

  \item{iso5}{differenciate trimming at 5 miRNA from rest}

  \item{iso3}{differenciate trimming at 3 miRNA from rest}

  \item{add}{differenciate additions miRNA from rest}

  \item{mism}{differenciate nt substitution miRNA from
  rest}

  \item{seed}{differenciate changes in 2-7 nt from rest}
}
\value{
count table
}
\description{
create count tables from isomirs
}

