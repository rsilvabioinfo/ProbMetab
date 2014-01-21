\name{incorporate.isotopes}
\alias{incorporate.isotopes}
\title{incorporate.isotopes}
\usage{
  incorporate.isotopes(ionAnnot, reactionM, comb = NULL,
    polarity = NULL, var = 2, samp = NULL, DB)
}
\arguments{
  \item{ionAnnot}{annotation list from get.annot function.}

  \item{reactionM}{compound's reaction matrix.}

  \item{comb}{1 for acquisition mode combination.}

  \item{polarity}{acquisition mode polarity.}

  \item{var}{1 to use standard mean/sd estimators to carbon
  number prediction, 2 for median/mad estimators.}

  \item{samp}{sample indexes, other than blanks, controls
  and QCs, according to xcms's phenoData.}

  \item{DB}{data.frame of compound information, with
  chemical formula.}
}
\value{
  matrix of candidate compound theoretical isotope
  patterns.
}
\description{
  Calculates the theoretical pattern of first 13C isotope
  for each candidate formula.
}

