\name{createJSONToCytoscape}
\alias{createJSONToCytoscape}
\title{createJSONToCytoscape}
\usage{
  createJSONToCytoscape(gr, node.label, node.form = NULL,
    edge.form = NULL, saveAsJSONFile = TRUE)
}
\arguments{
  \item{gr}{a graphNEL Object.}

  \item{node.label}{a vector of node labels (usually in the
  same node(gr) order).}

  \item{node.form}{data.frame file to set node parameters
  to cytoscape.js. First column should refers to node.names
  (or node.label) and will be ignored (as we use node.label
  to set node labels).  Every other column must be a
  logical one (i.e. TRUE/FALSE values per row), and it's
  name should be separeted by a period (.) (e.g.
  shape.triangle or color.#FFFF00) One need to mind only
  valid parameters will be correctly rendered by cytoscape
  web.  Valid parameres are: - width - color}

  \item{edge.form}{data.frame file to set edge parameters
  to cytoscape.js. First column must refers to edge.names
  (or edge.label) and it's value should be set as
  source~target (e.g. 819~821) (same pattern used for
  RCytoscape edgeNames() function) Every other column must
  be a logical one (i.e. TRUE/FALSE values per row), and
  it's name should be separeted by a period (.) (e.g.
  width.5 or color.#FFFF00) One need to mind: only valid
  parameters will be correctly rendered by cytoscape web.
  Valid parameres are: - shape (valid values: {rectangle,
  roundrectangle, ellipse, triangle, pentagon, hexagon,
  heptagon, octagon} ) - height - color}

  \item{saveAsJSONFile}{if TRUE (default) saves a file
  named "network.json" user's current R Working Directory.}
}
\value{
  createJSONToCytoscape returns a JSON-like string object.
}
\description{
  createJSONToCytoscape converts a graph-like file into a
  JSON-like file structure. It is mostly used to enable
  cytoscape.js render a graph.
}

