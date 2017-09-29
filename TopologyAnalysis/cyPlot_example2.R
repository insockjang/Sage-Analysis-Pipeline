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

# Here, we will use the cyPlot function to create a graphNEL object with nodes, edges, and their
# attributes from the corresponding node and edge data frames.

node.tbl <- data.frame(person = igraph::V(gD)$name, degree = igraph::V(gD)$degree, betweenness = igraph::V(gD)$betweenness)
edge.tbl <- dataSet.ext
colnames(edge.tbl) <- c("person_1", "person_2", "weight", "similarity")

# gD.cyt <- RCy3::cyPlot(node.tbl, edge.tbl)
# There were some issues with the original cyPlot function - it was assumed that there will be more than one
# data frame column of the character type or at least one numeric type column,
# and the lack of integer to numeric type conversion for edges, so I made changes and called a new function
# cyPlot_mod function
gD.cyt <- cyPlot_mod(node.tbl, edge.tbl)

# Next, we will create a new graph window in cytoscape
gDCW <- RCy3::CytoscapeWindow("Les Miserables using cyPlot", graph = gD.cyt, overwriteWindow = TRUE)

# We can display graph, with defaults color/size scheme
RCy3::displayGraph(gDCW)

# Now let's send/load node and edge attributes into Cytoscape
# Using the cyPlot function, all attributes are loaded into Cytoscape directly and we don't
# have to use additional functions to initialize or set node attributes

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
# RCy3::setLayoutProperties (gDCW, hlp[10], list (attraction_multiplier = 'similarity', nIterations = 1000))
# RCy3::layoutNetwork(gDCW, hlp[10])

# The layout looks different than before and does not change significalty based on whether we select gravity
# or attraction multiplier as its property

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
# Thus, this approach does not solve this issue

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