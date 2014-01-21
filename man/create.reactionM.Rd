\name{create.reactionM}
\alias{create.reactionM}
\title{create.reactionM}
\usage{
  create.reactionM(DB, molIon, ppm.tol)
}
\arguments{
  \item{DB}{dataframe containing the mandatory fields id,
  formula and reactions.}

  \item{molIon}{list of annotations provided by get.annot
  function.}

  \item{ppm.tol}{parts per million mass tolerance allowed
  in the mass search.}
}
\value{
  A matrix of reactions that each compound candidate,
  inside the mass window, can participate in the
  metabolism.
}
\description{
  This function matches a mass vector against a user
  provided database, inside an user provided mass tolerance
  window.
}

