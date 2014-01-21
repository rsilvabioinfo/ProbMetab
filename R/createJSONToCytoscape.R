#' createJSONToCytoscape
#' 
#' createJSONToCytoscape converts a graph-like file into a JSON-like file structure.
#' It is mostly used to enable cytoscape.js render a graph.
#' 
#' @param gr a graphNEL Object. 
#' @param node.label a vector of node labels (usually in the same node(gr) order).
#' 
#' @param node.form data.frame file to set node parameters to cytoscape.js. First column should refers to node.names (or node.label) 
#' and will be ignored (as we use node.label to set node labels).
#'    Every other column must be a logical one (i.e. TRUE/FALSE values per row), and it's name should be separeted by a period (.) (e.g. shape.triangle or color.#FFFF00)
#'    One need to mind only valid parameters will be correctly rendered by cytoscape web.
#'    Valid parameres are:
#'      - width
#'      - color
#' 
#' @param edge.form data.frame file to set edge parameters to cytoscape.js. First column must refers to edge.names (or edge.label) 
#'    and it's value should be set as source~target (e.g. 819~821) (same pattern used for RCytoscape edgeNames() function)
#'    Every other column must be a logical one (i.e. TRUE/FALSE values per row), and it's name should be separeted by a period (.) (e.g. width.5 or color.#FFFF00)
#'    One need to mind: only valid parameters will be correctly rendered by cytoscape web.
#'    Valid parameres are:
#'      - shape (valid values: {rectangle, roundrectangle, ellipse, triangle, pentagon, hexagon, heptagon, octagon} )
#'      - height
#'      - color
#'    
#' @param saveAsJSONFile if TRUE (default) saves a file named "network.json" user's current R Working Directory.
#' @return 
#'    createJSONToCytoscape returns a JSON-like string object.
#'
#' @export


createJSONToCytoscape = function(gr, node.label, node.form=NULL, edge.form=NULL, saveAsJSONFile=TRUE){
 # require(rjson)
  #### Aux functions
  format.form = function (edge.form){ # generic - works for any .form
    for(i in 2:dim(edge.form)[2]){#dim(edge.form)[2]
      name = strsplit(colnames(edge.form)[i], "\\.")[[1]][1];
      #def = NULL; if(name == "width") {def = 1}; if(name == "color") {def = "#ff1234"}
      content = strsplit(colnames(edge.form)[i], "\\.")[[1]][2];
      edge.form[which(edge.form[i] == TRUE),i] = content;
      edge.form[which(edge.form[i] == FALSE),i] = "#808080";
      colnames(edge.form)[i] <- name;
    }
    return(edge.form)
  }
  
  addFormatToJSON = function(edge.form, edge_){
    for(i in 2:dim(edge.form)[2]){
      name = strsplit(colnames(edge.form)[i], "\\.")[[1]][1];
      for(j in 1:dim(edge.form)[1]){
        edge_[[j]]$data[name] <- edge.form[j,i]
      }
    }
    return(edge_)
  }
  
  
  ## null --> node.form && edge.form : zero (use gr and node.label)
  ## null --> node.form && not null edge.form : one (use edge.form)
  ## not null --> node.form && null edge.form : two (use node.form)
  ## not null --> node.form && not null edge.form : three (use edge.form and node.form)
  t <- "error";
  if( (is.null(node.form))  && (is.null(edge.form)) ) {
    t <- "zero"
  }
  if( (is.null(node.form))  && (!is.null(edge.form)) ) {
    t <- "one"
  }
  if( (!is.null(node.form))  && (is.null(edge.form)) ) {
    t <- "two"
  }
  if( (!is.null(node.form))  && (!is.null(edge.form)) ) {
    t <- "three"
  }
  switch(t, 
         zero={
           print("zero")
           ### (use gr and node.label)
           node.label = gsub("\n", " ", node.label)
           node_ = lapply(node.label, function(x) { x})
           i<- 1; node_ = lapply(nodes(gr), function(x) { i <<- i+1; list(data=list(id=x, label=node.label[[i-1]])) } )
           names(node_) <- NULL
           
           edge_ = lapply(edgeNames(gr), function(x)  list(data=list(label=strsplit(as.character.default(x), "~")[[1]][2], source=strsplit(as.character.default(x), "~")[[1]][1], target=strsplit(as.character.default(x), "~")[[1]][2])))
           
           final_ = list(data=(list(nodes=node_, edges=edge_)))
           r = toJSON(final_)
           #return(r)
         },#end-of zero-case
         one={
           print("one")
           ## (use edge.form)
           
           node.label = gsub("\n", " ", node.label)
           node_ = lapply(node.label, function(x) { x})
           i<- 1; node_ = lapply(nodes(gr), function(x) { i <<- i+1; list(data=list(id=x, label=node.label[[i-1]])) } )
           names(node_) <- NULL
           
           edge.form = format.form(edge.form)
           edge_ = apply(edge.form,1, function(x) { list(data=list(label=x[[1]], source=strsplit(as.character.default(x[1]), "~")[[1]][1], target=strsplit(as.character.default(x[1]), "~")[[1]][2]    )   ) })#data=(lapply(x[-1],function(y) unlist(unlist(y))) )
           names(edge_) <- NULL
           edge_ = addFormatToJSON(edge.form, edge_)
           
           final_ = list(data=(list(nodes=node_, edges=edge_)))
           r = toJSON(final_)
           #return(r)
           
         },#end-of one-case
         two={
           print("two")
           ## (use node.form)
           node.form = format.form(node.form)
           node.label = gsub("\n", " ", node.label)
           node_ = lapply(node.label, function(x) { x})
           i<- 1; node_ = lapply(nodes(gr), function(x) { i <<- i+1; list(data=list(id=x, label=node.label[[i-1]])) } )
           names(node_) <- NULL
           node_ = addFormatToJSON(node.form, node_)
           
           edge_ = lapply(edgeNames(gr), function(x)  list(data=list(label=strsplit(as.character.default(x), "~")[[1]][2], source=strsplit(as.character.default(x), "~")[[1]][1], target=strsplit(as.character.default(x), "~")[[1]][2])))
           final_ = list(data=(list(nodes=node_, edges=edge_)))
           r = toJSON(final_)
           #return(r)
         },#end-of two-case
         three={
           print("three")
           ## (use edge.form and node.form)
           node.form = format.form(node.form)
           node.label = gsub("\n", " ", node.label)
           node_ = lapply(node.label, function(x) { x})
           i<- 1; node_ = lapply(nodes(gr), function(x) { i <<- i+1; list(data=list(id=x, label=node.label[[i-1]])) } )
           names(node_) <- NULL
           node_ = addFormatToJSON(node.form, node_)
           
           edge.form = format.form(edge.form)
           edge_ = apply(edge.form,1, function(x) { list(data=list(label=x[[1]], source=strsplit(as.character.default(x[1]), "~")[[1]][1], target=strsplit(as.character.default(x[1]), "~")[[1]][2]    )   ) })#data=(lapply(x[-1],function(y) unlist(unlist(y))) )
           names(edge_) <- NULL
           edge_ = addFormatToJSON(edge.form, edge_)
           
           final_ = list(data=(list(nodes=node_, edges=edge_)))
           r = toJSON(final_)
           #return(r)
         } #end-of three-case
  )#end-of switch
  if(saveAsJSONFile){
    write(paste("var json=JSON.parse('",r, "');", sep=""), file="network.json")
    print('File saved')
  }
  else {cat(r)}
}# end-of function
