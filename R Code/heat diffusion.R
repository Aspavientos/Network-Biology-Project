# Information ----
# Author: Diego Rodríguez
# Date of creation: 13/10/2023
# Last update: 16/10/2023

# Loading ----
## Load libraries ----
require(RCy3)
require(igraph)
require(clusterProfiler)
require(org.Hs.eg.db)
require(R.utils)
require(tidyverse)
require(sccore)

## Load files ----
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("..")

if(!file.exists("./PPI Network/9606.protein.links.v12.0.txt")){
  gunzip("./PPI Network/9606.protein.links.v12.0.txt.gz", remove=FALSE)
}

data <- read.table(gzfile("./PPI Network/9606.protein.links.v12.0.txt"), sep=" ", header=TRUE, stringsAsFactors=FALSE)

# Graph cleaning ----

# Cytoscape/igraph ----
## Connecting to Cytoscape ----

## igraph ----
graph <- graph_from_data_frame(data)

# Calculate heat diffusion ----
seed_prot = read.csv('./PPI Network/STRING network - endometriosis IDs.csv', header = FALSE, stringsAsFactors = FALSE)

seed_vec = seed_prot[,1]
names(seed_vec) = seed_vec

## sccore method ----
sccore_prop = propagate_labels(E(graph), rep(1, times = 13714404), seed_vec)

## Manually ----
# h
seed_ids = match(seed_prots$V1, as_ids(V(graph)))

# D
degs = degree(graph)

# A
adj_graph = as_adjacency_matrix(graph, sparse = FALSE)


# Diffusion
diffed = h*expm(-(degs-adj_graph)*0.1)

rm(degs, adj_graph)