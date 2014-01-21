\name{reac2cor}
\alias{reac2cor}
\title{reac2cor}
\usage{
  reac2cor(mw, classTable, opt = "cor", corths = 0.75,
    corprob = 0.8, pcorprob = 0.8)
}
\arguments{
  \item{mw}{two column of adjacency matrix indexes
  connecting compounds by reactions.}

  \item{classTable}{classification table, with intensities
  for repeated samples.}

  \item{opt}{correlation option, "cor" for correlation, and
  "pcor" for partial correlation.}

  \item{corths}{correlation intensity threshold.}

  \item{corprob}{probability that the correlation is
  considered significant.}

  \item{pcorprob}{probability that the partial correlation
  is considered significant.}
}
\value{
  A list of estimated correlations and reactions.
}
\description{
  Use the intensity of putative molecules in repeated
  samples to calculate correlations and partial correlation
  in a user defined threshold of false discovery rate for
  significance testing. After the correlation test the
  function also overlay significant correlations with all
  putative reactions between two masses.
}

