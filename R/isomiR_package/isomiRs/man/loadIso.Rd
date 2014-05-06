\name{loadIso}
\alias{loadIso}
\title{load data into IsomirDataSeq object
See Vignete to learn how to use this function}
\usage{
loadIso(files, design, cov = 1, header = F, skip = 1)
}
\arguments{
  \item{files}{all samples}

  \item{cov}{remove sequences that have relative abundance
  lower than this number}

  \item{design}{data frame containing groups for each
  sample}

  \item{header}{files contain headers}

  \item{skip}{skip first line when reading files}
}
\value{
\code{IsomirDataSeq} object
}
\description{
load data into IsomirDataSeq object See Vignete to learn
how to use this function
}

