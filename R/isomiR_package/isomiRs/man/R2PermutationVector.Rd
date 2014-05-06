\name{R2PermutationVector}
\alias{R2PermutationVector}
\title{Function that computes p-values when only individuals are considered to be permuted}
\usage{
R2PermutationVector(variables, group, validation, learn, test, tol, nperm)
}
\arguments{
  \item{variables}{matrix or data frame with explanatory
  variables}

  \item{group}{vector or factor with group memberships}

  \item{validation}{type of validation, either NULL or
  "learntest". Default NULL}

  \item{learn}{optional vector of indices for a learn-set.
  Only used when validation="learntest". Default NULL}

  \item{test}{optional vector of indices for a test-set.
  Only used when validation="learntest". Default NULL}

  \item{tol}{tolerance value based on maximum change of
  cumulative R-squared coefficient for each additional PLS
  component. Default tol=0.001}

  \item{nperm}{number of permutations to compute the PLD-DA
  p-value based on R2 magnitude. Default nperm=400}
}
\value{
PLS model
}
\description{
Function that computes p-values when only individuals are
considered to be permuted
}

