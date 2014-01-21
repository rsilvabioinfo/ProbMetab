\name{export.class.table}
\alias{export.class.table}
\title{export.class.table}
\usage{
  export.class.table(gibbsL = NULL, reactionM,
    molIon = NULL, probM = NULL, html = FALSE,
    filename = "test", burnIn = 3000, linkPattern = "kegg",
    m.test = "none", class1 = NULL, class2 = NULL,
    norm = FALSE, DB, prob = "count")
}
\arguments{
  \item{gibbsL}{a list of attributions and probabilities
  from gibbs.samp function.}

  \item{reactionM}{data.frame with compound annotation
  information.}

  \item{molIon}{non redundant ion annotation.}

  \item{probM}{optionally to gibbsL, a matrix of
  likelihoods.}

  \item{html}{logical, indicating whether a html file
  should be generated. This parameter uses the raw data to
  plot EICs and may be time consuming.}

  \item{filename}{html file name, the default is "test".}

  \item{burnIn}{how many samples of the gibbs sampler
  should be discarded.}

  \item{linkPattern}{which pattern should be linked to
  compound id, for now we have implemented "kegg",
  "pubchem" and "chebi" patterns.}

  \item{m.test}{statistical test to compare mean
  differences. This option is only available to single
  acquisition mode analysis, with options "t.test" and
  "anova".}

  \item{class1}{if the m.test is "t.test" first class to
  compare in the test, according with xcmsSet phenoData.}

  \item{class2}{if the m.test is "t.test" second class to
  compare in the test, according with xcmsSet phenoData.}

  \item{norm}{logical, if TRUE performs median
  normalization from (Anal. Chem. 2011, 83, 5864-5872).}

  \item{DB}{data.frame table used to search compounds, with
  the field name to be incorporated in the html table.}

  \item{prob}{how to calculate the probability to attribute
  a mass to a compound. Default is "count", which divide
  the number of times each identity was was attributed by
  the number of samples. Optionally the user could choose
  to use the mean of the probabilities of the identity,
  "mean".}
}
\value{
  A list with a matrix "classTable" with attributions and
  probabilities and indexes of selected masses from xcms
  peak table.
}
\description{
  Builds a matrix with the probability for all mass to
  candidate compounds assignments, by averaging the number
  of assignments obtained by the gibbs sampler algorithm or
  ordering the compound candidates with the likelihood
  matrix.
}

