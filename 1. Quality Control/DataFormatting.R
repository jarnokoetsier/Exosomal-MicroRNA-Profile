
#==============================================================================#

#   File:     DataFormatting.R
#   Author:   Jarno Koetsier
#   Date:     February 6, 2023

#==============================================================================#


###############################################################################

# 0. Preparation

###############################################################################

# Load packages
library(tidyverse)
library(readxl)

# Set working directory
setwd("E:/RTTproject/ExosomeAnalysis")

# Get data directory
dataDir <- "E:/RTTproject/OriginalData/Exosomes/Analysis_result/Summary/Final_Report/3_Summary/2_miRNAExpression"

# Get output directory
outputDir <- "E:/RTTproject/ExosomeAnalysis"


###############################################################################

# 1. Data collection and formatting

###############################################################################


# Load miRNA expression data
mirExpr <- read_excel(paste0(dataDir, "/Table5_All_Expressed_miRNA.xlsx"))

# Set column names
colnames(mirExpr) <- mirExpr[47,]

# Remove redudant rows
mirExpr <- mirExpr[48:nrow(mirExpr),]

# Filter for known miRNAs (gp1)
mirExpr <- mirExpr[mirExpr$group == "gp1",]

# Check if miRNA index is a unique identifier
length(unique(mirExpr$miRNA_Index)) == nrow(mirExpr)

# Get raw expression
exprMatrix <- as.matrix(mirExpr[,str_detect(colnames(mirExpr), "(raw)")])
exprMatrix <- apply(exprMatrix, 2, as.numeric)
rownames(exprMatrix) <- mirExpr$miRNA_Index
colnames(exprMatrix) <- str_remove_all(colnames(exprMatrix),"raw")
colnames(exprMatrix) <- str_remove_all(colnames(exprMatrix),"[()]")

# Save data
save(exprMatrix, file = paste0(outputDir,"exprMatrix_raw.RData"))

# Get feature data
featureInfo <- as.data.frame(mirExpr[,1:23])
save(featureInfo, file = paste0(outputDir,"featureInfo.RData"))

# Get sampleInfo
sampleInfo <- read.delim(paste0(outputDir,"sampleInfo.txt"))
all(colnames(exprMatrix) == sampleInfo$SampleID)
save(sampleInfo, file = paste(outputDir,"sampleInfo.RData"))





