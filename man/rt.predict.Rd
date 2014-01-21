\name{rt.predict}
\alias{rt.predict}
\title{rt.predict}
\usage{
  rt.predict(testSet, predSet, descData = NULL,
    descList = NULL, columnPar = NULL, voidTime = NULL)
}
\arguments{
  \item{testSet}{is a training set provided by user.}

  \item{predSet}{is a set of candidate compound to have
  their times predicted.}

  \item{descData}{is the user provided list of matrices of
  compound descriptors.}

  \item{descList}{is the optional user provided list of
  descriptors name calculated by rcdk package.}

  \item{columnPar}{is the optional user provided list of
  column parameters.}

  \item{voidTime}{column "dead" time.}
}
\value{
  A list of prediction parameters.
}
\description{
  Predict retention time from a given training set, and a
  matrix of compound descriptors.
}

