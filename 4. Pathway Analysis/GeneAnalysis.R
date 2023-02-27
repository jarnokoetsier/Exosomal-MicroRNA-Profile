#==============================================================================#

#   File:     PathwayAnalysis.R
#   Author:   Jarno Koetsier
#   Date:     February 6, 2023

#==============================================================================#


genes <- result_WP[1,"geneID"]

# "595/1026/207/23013/9794/3065/4790/2932/6774/3280/3714/4853/2033/2625/8819/
# 4609/83737/51107/6868/55294/3066/5295/3516/3091"
###############################################################################

# 0. Preparation

###############################################################################

# Load packages
library(tidyverse)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rWikiPathways)

# Set working directory
setwd("E:/RTTproject/ExosomeAnalysis/4. Pathway Analysis/")

# Get data directory
dataDir <- "E:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "E:/RTTproject/ExosomeAnalysis/4. Pathway Analysis/"


###############################################################################

# Data Preparation

###############################################################################

# Load mirRNA-target interaction data
load("MTI_filtered.RData")

# link gene to miRNA
gene2mir <- MTI_filtered[,c("Target Gene", "miRNA")]


# Get the relevent gene sets (Run only the relevent secton):

#==============================================================================#
# Early expressed miRNAs
#==============================================================================#

# Get early expressed miRNAs
miRNA_set <- read.table("E:/RTTproject/ExosomeAnalysis/2. Time Analysis/DownExpr.txt", col.names = "miRNA")


#==============================================================================#
# Lately expressed miRNAs
#==============================================================================#

# Get lately expressed miRNAs
miRNA_set <- read.table("E:/RTTproject/ExosomeAnalysis/2. Time Analysis/UpExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]

#==============================================================================#
# Ubiquitously expressed miRNAs
#==============================================================================#

# Get ubiquitously expressed miRNAs
miRNA_set <- read.table("E:/RTTproject/ExosomeAnalysis/2. Time Analysis/ConstantExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]

#==============================================================================#
# Differentially expressed miRNAs (RTT vs IC)
#==============================================================================#

# Get differentially expressed miRNAs (RTT vs IC
miRNA_set <- read.table("E:/RTTproject/ExosomeAnalysis/3. RTT vs IC/RTTvsIC.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]


###############################################################################

# WikiPathays Overrepresentation Analysis

###############################################################################

# Perform ORA
WP <- enricher(gene = miRNA_set$miRNA,
               universe = as.character(unique(gene2mir$miRNA)),
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               TERM2GENE = gene2mir)

result_WP <- WP@result
result_WP_copy <- WP@result[,c(1,2,5)]

# Number of permutations
nPerm <- 10000

# All miRNAs
all_mirnas <- unique(MTI_filtered$miRNA)

# Number of miRNAs to select during each permutation
n_miRNA <- length(miRNA_set$miRNA)