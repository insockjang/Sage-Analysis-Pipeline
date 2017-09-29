############################################################################################
############################################################################################
# Plotting networks in R - an example how to plot a network and 
# customize its appearance in Cytoscape directly from R using 
# the RCy3 package
############################################################################################
############################################################################################
# Clear workspace 
# rm(list = ls())
############################################################################################

# Read a data set. 
# Data format: dataframe with 3 variables; variables 1 & 2 correspond to interactions; variable 3 is weight of interaction
dataSet <- read.table("lesmis.txt", header = FALSE, sep = "\t")

# Create a graph. Use simplyfy to ensure that there are no duplicated edges or self loops
gD <- igraph::simplify(igraph::graph.data.frame(dataSet, directed=FALSE))

# Print number of nodes and edges
# igraph::vcount(gD)
# igraph::ecount(gD)

############################################################################################
# Calculate some node properties and node similarities that will be used to illustrate 
# different plotting abilities

# Calculate degree for all nodes
degAll <- igraph::degree(gD, v = igraph::V(gD), mode = "all")

# Calculate betweenness for all nodes
betAll <- igraph::betweenness(gD, v = igraph::V(gD), directed = FALSE) / (((igraph::vcount(gD) - 1) * (igraph::vcount(gD)-2)) / 2)
betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
rm(betAll)

#Calculate Dice similarities between all pairs of nodes
dsAll <- igraph::similarity.dice(gD, vids = igraph::V(gD), mode = "all")

############################################################################################
# Add new node and edge attributes based on the calculated node properties/similarities

gD <- igraph::set.vertex.attribute(gD, "degree", index = igraph::V(gD), value = degAll)
gD <- igraph::set.vertex.attribute(gD, "betweenness", index = igraph::V(gD), value = betAll.norm)

# Check the attributes
# summary(gD)

F1 <- function(x) {data.frame(V4 = dsAll[which(igraph::V(gD)$name == as.character(x$V1)), which(igraph::V(gD)$name == as.character(x$V2))])}
dataSet.ext <- plyr::ddply(dataSet, .variables=c("V1", "V2", "V3"), function(x) data.frame(F1(x)))

gD <- igraph::set.edge.attribute(gD, "weight", index = igraph::E(gD), value = 0)
gD <- igraph::set.edge.attribute(gD, "similarity", index = igraph::E(gD), value = 0)

# The order of interactions in dataSet.ext is not the same as it is in dataSet or as it is in the edge list
# and for that reason these values cannot be assigned directly

for (i in 1:nrow(dataSet.ext))
{
  igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$weight <- as.numeric(dataSet.ext$V3)
  igraph::E(gD)[as.character(dataSet.ext$V1) %--% as.character(dataSet.ext$V2)]$similarity <- as.numeric(dataSet.ext$V4)
}

# Check the attributes
# summary(gD)

rm(dataSet,dsAll, i, F1)

############################################################################################
# Now, let's do Cytoscape plots

# First, we need to transform our network from the igraph to graphnel format
gD.cyt <- igraph::as_graphnel(gD)

# Check if attributes have been passed
# graph::nodeData(gD.cyt, igraph::V(gD)$name, 'degree')
# graph::nodeData(gD.cyt, igraph::V(gD)$name, 'betweenness')
# graph::edgeData(gD.cyt, as.character(dataSet.ext$V1), as.character(dataSet.ext$V2), 'weight')
# graph::edgeData(gD.cyt, as.character(dataSet.ext$V1), as.character(dataSet.ext$V2), 'similarity')

# We have to create attributes for graphNEL
# We'll keep the same names as before
# In RCytoscape, this would ensure that the values of attributes are passed directly from igraph.
# However, this does not work with RCy3 right now (not sure if it is a bug or a feature has changed).
# Thus, we need to do send attributes to Cytoscape

gD.cyt <- RCy3::initNodeAttribute(gD.cyt, 'degree', 'numeric', 0) 
gD.cyt <- RCy3::initNodeAttribute(gD.cyt, 'betweenness', 'numeric', 0) 
gD.cyt <- RCy3::initEdgeAttribute (gD.cyt, "weight", 'integer', 0)
gD.cyt <- RCy3::initEdgeAttribute (gD.cyt, "similarity", 'numeric', 0)

# Next, we will create a new graph window in cytoscape
gDCW <- RCy3::CytoscapeWindow("Les Miserables", graph = gD.cyt, overwriteWindow = TRUE)

# We can display graph, with defaults color/size scheme
RCy3::displayGraph(gDCW)

# Now let's send/load node and edge attributes into Cytoscape

##########
# This should theoretically work, but there are some problems with attributes when networks
# are created from data frames (see https://github.com/tmuetze/Bioconductor_RCy3_the_new_RCytoscape/issues/25)
# I'll keep this code uncommented, but right now, it doesn't do anything

# setNodeAttributes should transfer the specified node attributes, for all nodes, the named node attribute 
# from the R graph (found in obj@graph) to Cytoscape. 
attribute.names <- RCy3::noa.names(gDCW@graph)

# Print list of attribute names to see if they are ok
# attribute.names 

# All nodes should already be in
RCy3::sendNodes(gDCW)

for (attribute.name in attribute.names){
  RCy3::setNodeAttributes(gDCW, attribute.name)
}

attribute.names <- RCy3::eda.names(gDCW@graph)

# All edges should already be in
RCy3::sendEdges(gDCW)

for (attribute.name in attribute.names){
  RCy3::setEdgeAttributes(gDCW, attribute.name)
}

RCy3::displayGraph(gDCW)

##########
# Thealternative, when we set attributes directly, works fine,
# so we will use it now (although, it seems kind of repetative)

RCy3::setNodeAttributesDirect(gDCW, 'degree', 'numeric', igraph::V(gD)$name, igraph::V(gD)$degree)
RCy3::setNodeAttributesDirect(gDCW, 'betweenness', 'numeric', igraph::V(gD)$name, igraph::V(gD)$betweenness)
RCy3::setEdgeAttributesDirect(gDCW, 'weight', 'integer', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(gD.cyt, as.character(dataSet.ext$V1), as.character(dataSet.ext$V2), 'weight'))
RCy3::setEdgeAttributesDirect(gDCW, 'similarity', 'numeric', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(gD.cyt, as.character(dataSet.ext$V1), as.character(dataSet.ext$V2), 'similarity'))

##########
# Now let's decide on a layout

# If you also want to choose a layout from R, a list  of available layouts can be accessed as follow:
cy <- RCy3::CytoscapeConnection()
hlp <-RCy3::getLayoutNames(cy)

# We'll select the "fruchterman-rheingold" layout. This layout is the layout number 10 
# To see properties for the given layout, use:
# RCy3::getLayoutPropertyNames(cy, hlp[10])
# We can choose any property we want and provide them as a list
RCy3::setLayoutProperties (gDCW, hlp[10], list (gravity_multiplier = 'similarity', nIterations = 1000))
RCy3::layoutNetwork(gDCW, hlp[10])

# I've noticed that if I change property to attraction_multiplier
# RCy3::setLayoutProperties (gDCW, hlp[10], list (attraction_multiplier = 'similarity', nIterations = 1000))
# RCy3::layoutNetwork(gDCW, hlp[10])
# And then go back to the original one
# RCy3::setLayoutProperties (gDCW, hlp[10], list (gravity_multiplier = 'similarity', nIterations = 1000))
# RCy3::layoutNetwork(gDCW, hlp[10])
# The layout won't go back to the original one. I am not sure if this is a bug or not

##########
# Finally, we can define rules for nodes:
RCy3::setNodeColorRule(gDCW, 'degree', c(min(degAll), mean(degAll), max(degAll)), c('#F5DEB3', '#FFA500', '#FF7F50', '#FF4500', '#FF0000'), mode = 'interpolate')
RCy3::setNodeSizeRule(gDCW, 'betweenness', c(min(betAll.norm), mean(betAll.norm), max(betAll.norm)), c(30, 45, 60, 80, 100), mode = 'interpolate')

# And edges:
RCy3::setEdgeLineWidthRule(gDCW, 'weight', dataSet.ext$V3, dataSet.ext$V3)
RCy3::setEdgeColorRule(gDCW, 'weight', c(min(as.numeric(dataSet.ext$V3)), mean(as.numeric(dataSet.ext$V3)), max(as.numeric(dataSet.ext$V3))), c('#FFFF00', '#00FFFF', '#00FF7F', '#228B22', '#006400'), mode='interpolate')

# While I get the "Successfully set rule" for both of the Edge rules, the view in the Cytoscape did not 
# change accordning the rules - setEdgeLineWidthRule command did not make any changes and the 
# setEdgeColorRule command made all edges white.
# One of the GitHub solved issues suggests to first set all rule-based functions and then the direct ones, but
# but it didn't work here (https://github.com/tmuetze/Bioconductor_RCy3_the_new_RCytoscape/issues/21 and
# https://github.com/tmuetze/Bioconductor_RCy3_the_new_RCytoscape/issues/20)

# We will define our own default color/size schema after we defined node and edge rules, due to
# possible issues when using rules
RCy3::setDefaultBackgroundColor(gDCW, '#D3D3D3')
RCy3::setDefaultEdgeColor(gDCW, '#CDC9C9')
RCy3::setDefaultEdgeLineWidth(gDCW, 4)
RCy3::setDefaultNodeBorderColor(gDCW, '#000000')
RCy3::setDefaultNodeBorderWidth(gDCW, 3)
RCy3::setDefaultNodeShape(gDCW, 'ellipse')
RCy3::setDefaultNodeColor(gDCW, '#87CEFA')
RCy3::setDefaultNodeSize(gDCW, 60)
RCy3::setDefaultNodeFontSize(gDCW, 20)
RCy3::setDefaultNodeLabelColor(gDCW, '#000000')