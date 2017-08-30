###################################################
### code chunk number 0: input data for g0
###################################################require(graph)
load("/gpfs/archive/RED/isjang/msigdb_v6.0_files_to_download_locally/MSigDBv6.Rdata")
geneSet <- MSigDB$H.ALL$genesets

JAK <- matrix(NA,nrow = length(geneSet),ncol = length(geneSet))
SIZE <- rep(NA,length(geneSet))

for(k1 in 1:length(geneSet)){
  SIZE[k1] <- length(geneSet[[k1]])
  if(k1 == length(geneSet)){
    next
  }else{
    for(k2 in (k1 + 1):length(geneSet)){
      # U<-min(length(geneSet[[k1]]),length(geneSet[[k2]]))
      # I<-length(intersect(geneSet[[k1]],geneSet[[k2]]))
      
      U<-length(union(geneSet[[k1]],geneSet[[k2]]))
      I<-length(intersect(geneSet[[k1]],geneSet[[k2]]))
      JAK[k1,k2]<-I/U
    }
  }
  print(k1)
}

#JAK[JAK == 0] <- NA
names(SIZE)<-rownames(JAK)<-colnames(JAK)<-gsub("HALLMARK_","",names(geneSet))


node.size <- SIZE
finalNet <- JAK

node.name <- names(node.size)

tab.HALLMARK <-read.delim("/gpfs/archive/RED/isjang/Project_CRISPR/JIS_CRISPR_OPM2/Analysis_normalizedNegControl/pathwayAnalysis/gsea_HALLMARK.txt", stringsAsFactors = F)
name.tab <- gsub("HALLMARK_","",tab.HALLMARK$geneset)
idx<-match(node.name, name.tab)
node.sig <-(tab.HALLMARK$fdr * sign(tab.HALLMARK$nes))[idx]

node.sig <-tab.HALLMARK$fdr[idx]

###################################################
### code chunk number 1: g0 (eval = FALSE)
###################################################
library(RCy3)

g <- new ('graphNEL', edgemode='undirected')
for(k in 1:length(node.name)){
  g <- graph::addNode (node.name[k], g)
}
cw <- CytoscapeWindow ('CRISPR_HALLMARK', graph=g, overwrite=TRUE)
displayGraph (cw)



###################################################
### code chunk number 2: g1 (eval = FALSE)
###################################################
## layoutNetwork (cw, layout.name='grid')
## 


###################################################
### code chunk number 3: g2 (eval = FALSE)
###################################################
g <- cw@graph   # created above, in the section  'A minimal example'
g <- initNodeAttribute (graph=g,  attribute.name='size',
                        attribute.type='numeric',
                        default.value=0.0)
for(k in 1:length(node.name)){
  nodeData (g, node.name[k], "size") <- node.size[k]
}

g <- initNodeAttribute (graph=g,  attribute.name='significance',
                        attribute.type='numeric',
                        default.value=0.0)

for(k in 1:length(node.name)){
  nodeData (g, node.name[k], "significance") <- node.sig[k]
}

cw <- setGraph (cw, g)
displayGraph (cw)    # cw's graph is sent to Cytoscape



###################################################
### code chunk number 4: defaults (eval = FALSE)
###################################################
setDefaultNodeShape (cw, 'ELLIPSE')
#setDefaultNodeColor (cw, rgb(1,1,0))
#setDefaultNodeSize  (cw, 40)
setDefaultNodeFontSize (cw, 15)


###################################################
### code chunk number 5: g3 (eval = FALSE)
###################################################
# getNodeShapes (cw)   # diamond, ellipse, trapezoid, triangle, etc.
# print (noa.names (getGraph (cw)))  # what data attributes are defined?
# print (noa (getGraph (cw), 'moleculeType'))
# attribute.values <- c ('kinase',  'TF',       'cytokine')
# node.shapes      <- c ('DIAMOND', 'TRIANGLE', 'RECTANGLE')
# setNodeShapeRule (cw, node.attribute.name='moleculeType',
#                   attribute.values, node.shapes)
# 
###################################################
### code chunk number 6: g4 (eval = FALSE)
###################################################
setNodeColorRule (cw, 'significance', seq(0,1,len = 10),
                  grey.colors(10, start = 0.9, end = 0.25, gamma = 2.2, alpha = NULL),
                  mode='interpolate')


###################################################
### code chunk number 8: g5 (eval = FALSE)
###################################################

setNodeSizeRule (cw, 'size', seq(min(node.size), max(node.size), length = 50), seq(min(node.size), max(node.size), length = 50)/3.5,
                 mode='interpolate')


#####
require(reshape2)
cutoff = 0.5
SDJ <- setNames(melt(finalNet), c('source', 'destination', 'Jakad'))
EDGE <- SDJ[which(SDJ$Jakad >= cutoff),]
###################################################
### code chunk number 9: g6 (eval = FALSE)
###################################################
g <- cw@graph
g <- initEdgeAttribute (graph=g,  attribute.name='Jakad',
                        attribute.type='numeric',
                        default.value=1)
g <- initEdgeAttribute (graph=g,  attribute.name='edgeType',
                        attribute.type='char',
                        default.value='unspecified')

for(k in 1:nrow(EDGE)){
  g <- graph::addEdge (as.character(EDGE$source[k]), as.character(EDGE$destination[k]), g)
}

for(k in 1:nrow(EDGE)){
  edgeData(g, as.character(EDGE$source[k]), as.character(EDGE$destination[k]), "Jakad") <- as.numeric(EDGE$Jakad[k])
  edgeData(g, as.character(EDGE$source[k]), as.character(EDGE$destination[k]), "edgeType") <- "indirected"
}


cw@graph <- g
displayGraph (cw)



setDefaultEdgeLineWidth(cw, 1)

edgeType.values = c ('indirected')

line.styles = c ('SOLID')
setEdgeLineStyleRule (cw, 'edgeType', edgeType.values,
                      line.styles)


setEdgeTargetArrowRule (cw, "edgeType", "indirected","NONE")

setEdgeLineWidthRule(cw, "Jakad", seq(0,1,length = 20),seq(0,20,length = 20))
