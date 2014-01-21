\name{create.reaction.edge.attributes}
\alias{create.reaction.edge.attributes}
\title{create.reaction.edge.attributes}
\usage{
  create.reaction.edge.attributes(classTable, graph, w,
    reactionM, DB, filename)
}
\arguments{
  \item{classTable}{classification table created by
  export.class.table function.}

  \item{graph}{graphNEL object, with node indexes
  corresponding to mass indexes in classTable.}

  \item{w}{matrix of compound connections.}

  \item{reactionM}{data.frame with compound annotation
  information.}

  \item{DB}{database with compound names associated to
  unique ids, used by export.class.table function.}

  \item{filename}{reaction attribute list file.}
}
\value{
  Writes a standard Cytoscape attribute list to the current
  working directory.
}
\description{
  This function writes a standard Cytoscape edge attribute
  list (http://www.cytoscape.org/) file. It takes the the
  possible mass connections, codified w matrix, retrieve
  the reactions where the compound is known to be present
  and associate then to a mass edge. It only works for
  KEGG, but, a specification of database will be available
  soon.
}

