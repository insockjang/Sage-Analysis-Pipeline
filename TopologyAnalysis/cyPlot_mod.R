cyPlot_mod <- function (node.df, edge.df) 
{
  edge.nodes <- unique(c(as.character(edge.df[, 1]), as.character(edge.df[, 
                                                                          2])))
  mydata <- new("graphNEL", edgemode = "directed", nodes = unique(c(as.character(node.df[, 
                                                                                         1]), edge.nodes)))
  node.df[, 1] <- as.character(node.df[, 1])
  edge.df[, 1:2] <- sapply(edge.df[, 1:2], as.character)
  node.class <- sapply(node.df, class)
  
  if (any(grep("factor", node.class))) {
    node.df[, grep("factor", node.class)] <- sapply(node.df[, 
                                                            grep("factor", node.class)], as.character)
  }
  
  if (any(grep("integer", node.class))) {
    node.df[, grep("integer", node.class)] <- sapply(node.df[, 
                                                             grep("integer", node.class)], as.numeric)
  }
  
  node.class <- sapply(node.df, class)
  edge.class <- sapply(edge.df, class)
  
  if (any(grep("factor", edge.class))) {
    edge.df[, grep("factor", edge.class)] <- sapply(edge.df[, 
                                                            grep("factor", edge.class)], as.character)
  }
  
  if (any(grep("integer", edge.class))) {
    edge.df[, grep("integer", edge.class)] <- sapply(edge.df[, 
                                                             grep("integer", edge.class)], as.numeric)
  }
  
  edge.class <- sapply(edge.df, class)
  
  
  if (length(grep("character", node.class)) > 1)
    for (i in 2:length(grep("character", node.class))) {
      mydata <- RCy3::initNodeAttribute(graph = mydata, attribute.name = names(node.class[grep("character", 
                                                                                               node.class)])[i], attribute.type = "char", default.value = "undefined")
      graph::nodeData(mydata, n = as.character(node.df[, 1]), attr = names(node.class[grep("character", 
                                                                                           node.class)])[i]) <- as.character(node.df[, grep("character", 
                                                                                                                                            node.class)[i]])
    }
  
  if (length(grep("numeric", node.class)) > 0)
    for (i in 1:length(grep("numeric", node.class))) {
      mydata <- RCy3::initNodeAttribute(graph = mydata, attribute.name = names(node.class[grep("numeric", 
                                                                                               node.class)])[i], attribute.type = "numeric", default.value = 0)
      graph::nodeData(mydata, n = as.character(node.df[, 1]), attr = names(node.class[grep("numeric", 
                                                                                           node.class)])[i]) <- as.numeric(node.df[, grep("numeric", 
                                                                                                                                          node.class)[i]])
    }
  
  mydata = graph::addEdge(as.vector(edge.df[, 1], mode = "character"), 
                          as.vector(edge.df[, 2], mode = "character"), mydata)
  
  if (length(grep("character", edge.class)) > 2)
    for (i in 3:length(grep("character", edge.class))) {
      mydata <- RCy3::initEdgeAttribute(graph = mydata, attribute.name = names(edge.df[, 
                                                                                       grep("character", edge.class)])[i], attribute.type = "char", 
                                        default.value = "undefined")
      graph::edgeData(mydata, as.vector(edge.df[, 1], mode = "character"), 
                      as.vector(edge.df[, 2], mode = "character"), attr = names(edge.df[, 
                                                                                        grep("character", edge.class)])[i]) <- as.character(edge.df[, 
                                                                                                                                                    grep("character", edge.class)[i]])
    }
  
  if (length(grep("numeric", edge.class)) > 0)
    for (i in 1:length(grep("numeric", edge.class))) {
      mydata <- RCy3::initEdgeAttribute(mydata, attribute.name = names(edge.class[grep("numeric", 
                                                                                       edge.class)])[i], attribute.type = "numeric", default.value = 0)
      graph::edgeData(mydata, as.vector(edge.df[, 1], mode = "character"), 
                      as.vector(edge.df[, 2], mode = "character"), attr = names(edge.class[grep("numeric", 
                                                                                                edge.class)])[i]) <- as.numeric(edge.df[, grep("numeric", 
                                                                                                                                               edge.class)[i]])
    }
  return(mydata)
}