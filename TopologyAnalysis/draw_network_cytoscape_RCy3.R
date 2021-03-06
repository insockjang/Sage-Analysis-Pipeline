drawEnrichMapCytoscape<-function(DB,subDB,NODE,size,stats,similarity.cutoff = 0.5){
  # NODE <- reactome geneset names 
  # size <- reactome geneset size 
  # stats <- reactome geneset's gsea significance score : -log10(fdr) * sign(nes)
  # drawEnrichMapCytoscape(mat,NODE,size,nes, similarity.cutoff = 0.5)
  
  if(toupper(DB) == "GRAPHITE"){
    githubURL <- paste0("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/TopologyAnalysis/graphite.",toupper(subDB),"_similarity.Rdata")
    load(url(githubURL))
  }
  
  if(toupper(DB) == "MSIGDB"){
    githubURL <- paste0("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/TopologyAnalysis/msigdb.",toupper(subDB),".ALL_similarity.Rdata")
    load(url(githubURL))
  }
  
  if(toupper(DB) == "CPDB"){
    githubURL <- paste0("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/TopologyAnalysis/cpdb_similarity.Rdata")
    load(url(githubURL))
  }
  
  require(RCy3)
  require(reshape2)
  # NODE<-nodes(rEG)
  
  g <- new ("graphNEL",edgemode="undirected")
  for(k in 1:length(NODE)){
    g <- graph::addNode (NODE[k],g)
  }
  
  cw <- CytoscapeWindow ("vignette",graph=g, overwrite=TRUE)
  displayGraph (cw)
  
  
  g <- cw@graph
  
  g <- initNodeAttribute (graph=g, "size", "numeric", 0.0)
  g <- initNodeAttribute (graph=g, "stats", "numeric", 0.0)
  
  for(k in 1:length(NODE)){
    nodeData (g,NODE[k],"size")<- size[k]
    nodeData (g,NODE[k],"stats")<- stats[k]
  }
  
  cw <- setGraph (cw, g)
  displayGraph (cw) 
  
  setDefaultNodeShape (cw,"ELLIPSE")
  setDefaultNodeColor (cw,"#AAFF88")
  
  setDefaultNodeSize  (cw, 80)
  setDefaultNodeFontSize (cw, 10)
  
  # getNodeShapes (cw)   # diamond, ellipse, trapezoid, triangle, etc.
  print (noa.names (getGraph (cw)))  # what data attributes are defined?
  print (noa (getGraph (cw),"size"))
  
  control.point <- c(-4,0,4)
  setNodeColorRule (cw, "stats", control.point,c ("#0000FF","#FFFFFF","#FF0000"),mode="interpolate")
  
  setNodeSizeRule (cw, 'label', NODE,  log2(size) ^1.75, default.size= 10, mode='lookup')
  
  displayGraph(cw)
  
  # setEdgeAttributes (cw, "edgeWidth")
  setDefaultEdgeLineWidth  (cw, 0.1)
  setDefaultEdgeColor(cw,"#d3d3d3")
  #g <- initEdgeAttribute (graph=g,  attribute.name="edgeName",attribute.type="char",default.value="undefined")
  g <- initEdgeAttribute (graph=g,  attribute.name="edgeWidth",attribute.type="numeric",default.value=0.0)
  
  mat <- sim.combine
  mat[lower.tri(mat, diag = T)] <- NA
  
  mat <- mat[NODE,NODE]
  A <- setNames(melt(mat), c('from', 'to', 'score'))
  
  id.na <- which(is.na(A$score))
  AA <- A[-id.na,]
  AAA <- AA[which(AA$score !=0),]
  tab<- AAA[which(AAA$score >= similarity.cutoff),]
  
  
  
  for(k in 1:nrow(tab)){
    g <- graph::addEdge (as.character(tab$from[k]), as.character(tab$to[k]),g)
  }
  
  # edgeName <- paste(as.character(tab$from),as.character(tab$to),sep = " (undefined) ")
  for(k in 1:nrow(tab)){
    edgeData (g,as.character(tab$from[k]), as.character(tab$to[k]),"edgeWidth") <- tab$score[k]
  }
  
  cw@graph <- g
  
  displayGraph (cw)
  
  layoutNetwork(cw, "cose")
  return(cw)
}


drawPathwayCytoscapeMSigDB <- function(pathway.msigdb, reference){
  
  # pathway should be from MSigDB reactome database 
  # curated by In Sock Jang, edgeType = "directed" or "undirected"
  # find mapping database in the following
  
  
  githubURL <- "https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/TopologyAnalysis/MSigDB_GRAPHITE.reactome_mapping.txt"
  mapTable <- read.delim(url(githubURL), header = F, stringsAsFactors = F)
  id <- which(mapTable$V1 %in% pathway.msigdb)
  
  githubURL1 <- "https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/Graphite.Rdata"
  load(url(githubURL1))
  pathway <- pathwayGraph(GRAPHITE$REACTOME[[mapTable$V2[id]]])
  
  
  require(RCy3)
  
  NODE<-nodes(pathway)
  
  id<-match(names(reference), NODE)
  statistics <- rep(0, length(NODE))
  names(statistics) <- NODE
  statistics[id[which(!is.na(id))]] <- reference[which(!is.na(id))]
  
  pathway <- initEdgeAttribute (pathway, "weight", "numeric", 1)
  pathway <- initEdgeAttribute (pathway, "edgeType", "char", "undefined")
  
  cw <- CytoscapeWindow ("vignette",graph=pathway, overwrite=TRUE)
  displayGraph (cw)
  
  
  pathway <- cw@graph
  
  pathway <- initNodeAttribute (pathway, "stats", "numeric", 0.0)
  
  for(k in 1:length(NODE)){
    nodeData (pathway,NODE[k],"stats")<- statistics[k]
  }
  
  cw <- setGraph (cw, pathway)
  displayGraph (cw) 
  
  setDefaultNodeShape (cw,"ELLIPSE")
  setDefaultNodeColor (cw,"#AAFF88")
  
  setDefaultNodeSize  (cw, 40)
  setDefaultNodeFontSize (cw, 10)
  
  setDefaultNodeBorderWidth (cw, 1)
  setDefaultNodeBorderColor (cw, "#C0C0C0") # white borders
  
  # getNodeShapes (cw)   # diamond, ellipse, trapezoid, triangle, etc.
  print (noa.names (getGraph (cw)))  # what data attributes are defined?
  print (noa (getGraph (cw),"stats"))
  
  control.point <- c(-0.4,0,0.4)
  setNodeColorRule (cw, "stats", control.point,c ("#0000FF","#FFFFFF","#FF0000"),mode="interpolate")
  
  #setNodeSizeRule (cw, 'label', NODE,  log2(size) ^1.75, default.size= 10, mode='lookup')
  
  displayGraph(cw)
  
  # setEdgeAttributes (cw, "edgeWidth")
  setDefaultEdgeLineWidth  (cw, 0.5)
  setDefaultEdgeColor(cw,"#d3d3d3")
  
  setDefaultEdgeSourceArrowColor(cw,"#d3d3d3")
  setDefaultEdgeTargetArrowColor(cw,"#d3d3d3")
  edgeType.values = names(table(eda(pathway,"edgeType")))
  
  line.styles = rep("SOLID",length(edgeType.values))
  line.styles[grep("undirected",edgeType.values)]<- "LONG_DASH"
  setEdgeLineStyleRule (cw,"edgeType", edgeType.values, line.styles)
  
  
  arrow.styles = rep("ARROW", length(edgeType.values))
  arrow.styles[grep("undirected",edgeType.values)]<- "NONE"
  setEdgeTargetArrowRule (cw,"edgeType",edgeType.values, arrow.styles)
  # setEdgeSourceArrowRule (cw,"edgeType",edgeType.values, arrow.styles)
  
  
  layoutNetwork(cw, "cose")
  return(cw)
}




drawPathwayCytoscapeAdvance<-function(pathway, reference){
  
  # pathway should be from adjusted reactome database in advanced way : ACTIVATION and INHIBITION on directed edges
  # curated by In Sock Jang, edgeType = "directed" or "undirected"
  # find database in the following folder
  # /gpfs/archive/RED/isjang/PathwayDB/directAdvancedReactome.Rdata
  require(RCy3)
  
  NODE<-nodes(pathway)
  
  id<-match(NODE, names(reference))
  statistics <- reference[id]
  
  
  pathway <- initEdgeAttribute (pathway, "weight", "numeric", 1)
  pathway <- initEdgeAttribute (pathway, "edgeType", "char", "undefined")
  
  cw <- CytoscapeWindow ("vignette",graph=pathway, overwrite=TRUE)
  displayGraph (cw)
  
  
  pathway <- cw@graph
  
  pathway <- initNodeAttribute (pathway, "stats", "numeric", 0.0)
  
  for(k in 1:length(NODE)){
    nodeData (pathway,NODE[k],"stats")<- statistics[k]
  }
  
  cw <- setGraph (cw, pathway)
  displayGraph (cw) 
  
  setDefaultNodeShape (cw,"ELLIPSE")
  setDefaultNodeColor (cw,"#AAFF88")
  
  setDefaultNodeSize  (cw, 40)
  setDefaultNodeFontSize (cw, 10)
  
  # getNodeShapes (cw)   # diamond, ellipse, trapezoid, triangle, etc.
  print (noa.names (getGraph (cw)))  # what data attributes are defined?
  print (noa (getGraph (cw),"stats"))
  
  control.point <- c(-0.4,0,0.4)
  setNodeColorRule (cw, "stats", control.point,c ("#0000FF","#FFFFFF","#FF0000"),mode="interpolate")
  
  #setNodeSizeRule (cw, 'label', NODE,  log2(size) ^1.75, default.size= 10, mode='lookup')
  
  displayGraph(cw)
  
  # setEdgeAttributes (cw, "edgeWidth")
  setDefaultEdgeLineWidth  (cw, 0.5)
  setDefaultEdgeColor(cw,"#d3d3d3")
  
  setDefaultEdgeSourceArrowColor(cw,"#d3d3d3")
  setDefaultEdgeTargetArrowColor(cw,"#d3d3d3")
  edgeType.values = names(table(eda(pathway,"edgeType")))
  
  line.styles = rep("SOLID",length(edgeType.values))
  line.styles[grep("undirected",edgeType.values)]<- "LONG_DASH"
  setEdgeLineStyleRule (cw,"edgeType", edgeType.values, line.styles)
  
  
  arrow.styles = rep("ARROW", length(edgeType.values))
  arrow.styles[grep("undirected",edgeType.values)]<- "NONE"
  setEdgeTargetArrowRule (cw,"edgeType",edgeType.values, arrow.styles)
  
  # getLineStyles(cw)
  # getArrowShapes(cw)
  
  # edge.name <- cy2.edge.names (cw@graph)
  
  # setEdgeLineWidthRule(cw,"name",edgeName,tab$score)
  # id<-match(edge.name,edgeName)
  # duplicated(id)
  # setEdgeLineWidthDirect(cw,"edgeType",tab$score[1:188])  
  # 
  # 
  # 
  # print (noa.names (getGraph (cw)))  # what data attributes are defined?
  # print (noa (getGraph (cw),"size"))
  # length(cy2.edge.names (cw@graph))
  # eda.names(pathway)
  # eda(pathway)
  # getEdgeAttributeNames(cw)
  # edgeData(cw@graph)
  # 
  layoutNetwork(cw, "cose")
  return(cw)
}


drawCytoscapeNetworkPropagation<-function(mat, F.tab, cutoff.npzscore, weight.cutoff = 0.5){
  # input:
  #       mat : adjacent matrix pruned by network propagation 
  #       seed : basic given seeds for network propagation such as "CRBN", "GSPT1", "IKZF1", "IKZF3", etc
  # copyrighted by In Sock Jang (isjang@celgene.com)
  
  # F.tab has to have f.final, f.init, zscore columns : data frame or matrix class
  if(nrow(F.tab) != nrow(mat)){
    error("Please check if your network propagation adjacent matrix size is equal to the outcome table")
    break;
  }
  
  seeds <- colnames(mat)[which(F.tab$f.init == 1)]
  
  adds <- setdiff(colnames(mat)[which(F.tab$zscore > cutoff.npzscore)],seeds)
  
  allnodes <- union(seeds, adds)
  
  require(RCy3)
  require(reshape2)
  
  # normalize
  temp.mat <- mat[allnodes,allnodes]/1000
  
  
  temp.mat[lower.tri(temp.mat, diag = T)] <- NA
  
  A <- setNames(melt(temp.mat), c('from', 'to', 'score'))
  
  id.remove <- which(is.na(A$score) | A$score < weight.cutoff)
  tab<- A[-id.remove,]
  
  NODE <- union(tab$from, tab$to)
  
  g <- new ("graphNEL",edgemode="undirected")
  for(k in 1:length(NODE)){
    g <- graph::addNode (NODE[k],g)
  }
  
  
  for(k in 1:nrow(tab)){
    g <- graph::addEdge (as.character(tab$from[k]), as.character(tab$to[k]),g)
  }
  
  cw <- CytoscapeWindow ("vignette",graph=g, overwrite=TRUE)
  
  displayGraph (cw)
  
  g <- cw@graph
  g <- initEdgeAttribute (graph=g,  attribute.name="edgeWidth",attribute.type="numeric",default.value=0.0)
  g <- initEdgeAttribute (g,  "edgeType",attribute.type="char","undefined")
  for(k in 1:nrow(tab)){
    edgeData (g,as.character(tab$from[k]), as.character(tab$to[k]),"edgeWidth") <- tab$score[k]
  }
  
  g <- initNodeAttribute (g, "seeds", "char", "neighbor")
  
  id1 <- match(seeds,NODE)
  nodeData (g,NODE[id1[!is.na(id1)]],"seeds")<- "seed"
  
  cw <- setGraph (cw, g)
  displayGraph (cw) 
   
  setDefaultNodeShape (cw,"ELLIPSE")
  setDefaultNodeColor (cw,"#AAFF88")
  
  setDefaultNodeSize  (cw, 40)
  setDefaultNodeFontSize (cw, 10)
  
  # getNodeShapes (cw)   # diamond, ellipse, trapezoid, triangle, etc.
  # print (noa.names (getGraph (cw)))  # what data attributes are defined?
  # print (noa (getGraph (cw),"seeds"))
  
  data.values <- c("seed", "neighbor")
  node.colors <- c("#FFA500", "#7FFF00") 
  setNodeColorRule (cw, node.attribute.name="seeds", data.values, node.colors, mode='lookup', default.color='#AAFF88')
  
  shapes <- c("ROUND_RECTANGLE", "ELLIPSE")
  setNodeShapeRule (cw, node.attribute.name="seeds", data.values, shapes)
  
  arrows <- c("NONE")
  edgeType.values <- c("undefined")
  setEdgeTargetArrowRule (cw, 'edgeType', edgeType.values, arrows)
  
  setDefaultEdgeLineWidth  (cw, 0.5)
  setDefaultEdgeColor(cw,"#d3d3d3")
  
  displayGraph(cw)
  layoutNetwork(cw, "cose")
  return(cw)
}
