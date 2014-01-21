\name{create.pathway.node.attributes}
\alias{create.pathway.node.attributes}
\title{create.pathway.node.attributes}
\usage{
  create.pathway.node.attributes(classTable, graph, DB,
    filename1, filename2 = NULL, organismId = NULL)
}
\arguments{
  \item{classTable}{classification table created by
  export.class.table function.}

  \item{graph}{graphNEL object, with node indexes
  corresponding to mass indexes in classTable.}

  \item{DB}{database with compound names associated to
  unique ids, used by export.class.table function.}

  \item{filename1}{filename to attribute pathway list
  file.}

  \item{filename2}{optional filename to attribute pathway
  list discriminating compound/pathway associations.}

  \item{organismId}{KEGG organism id
  (http://www.kegg.jp/kegg/catalog/org_list.html) to filter
  possibibly pathwyas for known pathways for that
  organism.}
}
\value{
  writes a standard Cytoscape attribute list to current
  working directory. Also creates a matrix containing
  putative compound counting to each pathway.
}
\description{
  This function writes a standard Cytoscape node attribute
  list (http://www.cytoscape.org/) file. It takes the
  compound codes and retrieve all the known pathways, where
  the compound is known to be present. It only works for
  KEGG, but, a specification of database will be available
  soon.
}

