\name{gibbs.samp}
\alias{gibbs.samp}
\title{gibbs.samp}
\usage{
  gibbs.samp(x, y, N, w, p)
}
\arguments{
  \item{x}{a vector of masses (unique from mass/retention
  time pairs).}

  \item{y}{a vector of candidate compounds for each mass.}

  \item{N}{number of iterations to sample.}

  \item{w}{matrix of compound connections.}

  \item{p}{matrix of likelihood probabilities.}
}
\value{
  A list of matrices with attribution indexes and
  probabilities.
}
\description{
  Call a C++ function that runs a Gibbs Sampler algorithm
  to sample from the distribution of metabolite to compound
  attribution with the previous assumption that the
  connected combination of attributions makes more sense.
}

