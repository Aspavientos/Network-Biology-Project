# Information ----
# Author: Diego Rodríguez
# Date of creation: 13/10/2023
# 

# Loading ----
# Load libraries
require(regnet)
require(RCy3)
require(igraph)
require(clusterProfiler)
require(org.Hs.eg.db)
require(R.utils)

# Load files
if(!file.exists("./PPI Network/9606.protein.links.v12.0.txt")){
  gunzip("./PPI Network/9606.protein.links.v12.0.txt.gz", remove=FALSE)
}

data <- read.table(gzfile("./PPI Network/9606.protein.links.v12.0.txt"), sep=" ", header=TRUE, stringsAsFactors=FALSE)

# Graph cleaning ----


# Cytoscape/igraph ----
# Connecting to Cytoscape

# igraph
graph <- graph_from_data_frame(data)

# Calculate heat diffusion ----
seed_nodes = read.csv('./PPI Network/STRING network - endometriosis IDs.csv', header = FALSE, stringsAsFactors = FALSE)

diffusion <- seeded.diffusion(graph, seed_nodes, steps = 10)
