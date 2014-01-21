\name{get.kgml.positions.kegg}
\alias{get.kgml.positions.kegg}
\title{get.kgml.positions.kegg}
\usage{
  get.kgml.positions.kegg(path)
}
\arguments{
  \item{path}{KEGG's pathway code.}
}
\value{
  a list with an adjacency matrix for pathway nodes, and a
  matrix with node positions.
}
\description{
  Gets KEGG's pathway code, download the kgml file, and
  retrieves the pathway layout.
}

