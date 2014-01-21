\name{get.Mzmatch.annot}
\alias{get.Mzmatch.annot}
\title{get.Mzmatch.annot}
\usage{
  get.Mzmatch.annot(filename, onlyBP = TRUE)
}
\arguments{
  \item{filename}{PeakML file containing ion annotation.}

  \item{onlyBP}{logical, if TRUE retrieves only PeakML "bp"
  relationship, if FALSE also retrieves "potential bp"
  relationship.}
}
\value{
  A list with a matrix of possible molecular ions with a
  trace of their annotation, and the used xsAnnotate
  object.
}
\description{
  This function extracts annotation from mzMatch PeakML
  file, generating a matrix of non-redundant putative
  molecular ions.
}

