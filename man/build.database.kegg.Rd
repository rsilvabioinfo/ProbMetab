\name{build.database.kegg}
\alias{build.database.kegg}
\title{build.database.kegg}
\usage{
  build.database.kegg(orgID = NULL)
}
\arguments{
  \item{orgID}{KEGG's organism id. If NULL the function
  recovers information from all database compounds.}
}
\value{
  A data.frame with unique id, name, formula and reactions
  as needed by ProbMetab.
}
\description{
  Get KEGG compound information needed for exact mass
  searching and modeling.
}

