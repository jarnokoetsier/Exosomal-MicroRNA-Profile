# Load packages
library(tidyverse)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rWikiPathways)
library(RCy3)

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Set working directory
setwd("D:/RTTproject/ExosomeAnalysis/4. Pathway Analysis/")

# Get data directory
dataDir <- "D:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "D:/RTTproject/ExosomeAnalysis/4. Pathway Analysis/"

# Load feature info
load(paste0(dataDir,"featureInfo.RData"))

# Put miRNA name in correct format
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")

# Remove duplicated features
featureInfo <- featureInfo[!duplicated(featureInfo$CompleteName),]
cluster_chr14 <- featureInfo[featureInfo$genomeID == "chr14",]
cluster_chr14$start <- as.numeric(cluster_chr14$start)
cluster_chr14 <- cluster_chr14[(cluster_chr14$start < 101100000) &
                                 (cluster_chr14$start > 100800000),]
cluster_chr14 <- cluster_chr14[cluster_chr14$strand == "+",]


# Load the Serotonin pathway from WikiPathways

#   Edges
pathway <- read.csv("SerotoninPathway.csv")
source <- str_remove(str_remove(pathway$name,"\\(.*"), " ")
target <- str_remove(str_remove(pathway$name,".*\\)"), " ")

#   Nodes
nodes_serotonin <- read.csv("SerotoninPathway_nodes.csv")
nodes <- nodes_serotonin[,c("name", "Type")]

# Load mirRNA-target interaction data

# Edges
load("MTI_filtered.RData")
miRNA_set <- read.table("D:/RTTproject/ExosomeAnalysis/3. RTT vs IC/RTTvsIC.txt", col.names = "miRNA")
MTI_filtered <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]
MTI_filtered <- MTI_filtered[MTI_filtered$`Target Gene` %in% unique(c(source,target),),]
MTI_fil <- unique(MTI_filtered[,c("miRNA", "Target Gene")])

# Nodes
miRDF <- data.frame(name = MTI_fil$miRNA,
                    Type = ifelse(MTI_fil$miRNA %in% cluster_chr14$miR_name, "miRNA (C14MC)", "miRNA (other)"))



# Combine miRNA-target info with the Serotonin pathway

# Edges
edges <- data.frame(source = c(source,MTI_fil$miRNA),
                    target = c(target,MTI_fil$`Target Gene`))

# Nodes
nodes <- rbind.data.frame(nodes, miRDF)
colnames(nodes) <- c("id", "type")

# Create network in cytoscape
createNetworkFromDataFrames(
  nodes = nodes,
  edges = edges,
  title = "From dataframe",
  collection = "My Dataframe Network Collection",
)

