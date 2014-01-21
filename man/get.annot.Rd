\name{get.annot}
\alias{get.annot}
\title{get.annot}
\usage{
  get.annot(xsAnnotate, polarity = "positive",
    allowMiss = FALSE, xset = NULL, toexclude = NULL,
    minsamp = 0.6, minint = 5000)
}
\arguments{
  \item{xsAnnotate}{CAMERA's annotation object.}

  \item{polarity}{acquisition mode of mass spectrometer.}

  \item{allowMiss}{logical, optionally retrieves peaks with
  no evidence of adduct or isotope and annotate them as
  single charged molecules [M+/-H].}

  \item{toexclude}{samples to be excluded of peak counting
  to non-annotated peak selection.}

  \item{xset}{xcmsSet xcms object after missing data
  replacement, to retrieve SNR to isotopic peaks.}

  \item{minsamp}{minimum number of samples in which an ion
  should be present to be selected.}

  \item{minint}{minimum mean intensity that a ion should
  present to be selected.}
}
\value{
  A list with a matrix of possible molecular ions with a
  trace of their annotation, and the used xsAnnotate
  object.
}
\description{
  This function extracts annotation from CAMERA object,
  generating a matrix of non-redundant putative molecular
  ions.
}

