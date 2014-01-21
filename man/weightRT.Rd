\name{weightRT}
\alias{weightRT}
\title{weightRT}
\usage{
  weightRT(rtObj, reactionM, userCuttoff = 0.95,
    rtWeigth = 0.1, plot = FALSE)
}
\arguments{
  \item{rtObj}{is a list of fitted values predicted by
  rt.predict.}

  \item{reactionM}{is the KEGG search list.}

  \item{userCuttoff}{is the user set expectation of error
  acceptance.}

  \item{rtWeigth}{is penalty to errors that fall bellow the
  threshold given b user.}

  \item{plot}{logical, wheter or not do plot the output of
  the function.}
}
\value{
  A matrix wrt of likelihood.
}
\description{
  Builds a c (number of compounds) by m (number of masses)
  matrix of weights
}

