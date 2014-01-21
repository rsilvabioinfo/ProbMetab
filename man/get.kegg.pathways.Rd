\name{get.kegg.pathways}
\alias{get.kegg.pathways}
\title{get.kegg.pathways}
\usage{
  get.kegg.pathways(cpds, numCpds)
}
\arguments{
  \item{cpds}{KEGG compound ids, including "cpd:".}

  \item{numCpds}{minimum number of compounds as nodes of
  the pathway.}
}
\value{
  list of the missing compounds of kegg pathways, and a
  table with map links and statistics.
}
\description{
  Get pathway map links to all kegg pathways that have more
  than numCpds as nodes.
}

