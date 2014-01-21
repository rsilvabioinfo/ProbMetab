\name{sbml2table}
\alias{sbml2table}
\title{sbml2table}
\usage{
  sbml2table(file)
}
\arguments{
  \item{file}{SBML file.}
}
\value{
  A database in the format required for ProbMetab, with
  mandatory fields id, formula and reactions.
}
\description{
  Retrieves compound and its reactions from SBML models.
}

