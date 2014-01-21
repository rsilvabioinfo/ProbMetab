\name{comp.results}
\alias{comp.results}
\title{comp.results}
\usage{
  comp.results(reactionM, w, ansLik, ansConn)
}
\arguments{
  \item{reactionM}{matrix with reactions of each candidate
  compound.}

  \item{w}{matrix of compound connections.}

  \item{ansLik}{a list of mass to compound assignment,
  based only on a likelihood.}

  \item{ansConn}{list of mass to compound assignment, with
  compound's connections contribution.}
}
\value{
  a list with compound classification classes, table index
  of classes and a matrix of intensities of selected
  compounds.
}
\description{
  Compare two classification tables, given by
  export.class.table, and reports the difference between
  two different models.
}

