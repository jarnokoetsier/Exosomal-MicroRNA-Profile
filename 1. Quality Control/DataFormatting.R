# load packages
library(tidyverse)
library(readxl)

# Set working directory
setwd("E:/RTTproject/ExosomeAnalysis")

# Load miRNA expression data
dataDir <- "E:/RTTproject/OriginalData/Exosomes/Analysis_result/Summary/Final_Report/3_Summary/2_miRNAExpression"
mirExpr <- read_excel(paste0(dataDir, "/Table5_All_Expressed_miRNA.xlsx"))

# Set column names
colnames(mirExpr) <- mirExpr[47,]

# Remove redudant rows
mirExpr <- mirExpr[48:nrow(mirExpr),]

# Filter for known miRNAs
mirExpr <- mirExpr[mirExpr$group == "gp1",]

# Check if miRNA index is good identifier
length(unique(mirExpr$miRNA_Index)) == nrow(mirExpr)

# Get raw expression
exprMatrix <- as.matrix(mirExpr[,str_detect(colnames(mirExpr), "(raw)")])
exprMatrix <- apply(exprMatrix, 2, as.numeric)
rownames(exprMatrix) <- mirExpr$miRNA_Index
colnames(exprMatrix) <- str_remove_all(colnames(exprMatrix),"raw")
colnames(exprMatrix) <- str_remove_all(colnames(exprMatrix),"[()]")

# Save data
save(exprMatrix, file = "exprMatrix_raw.RData")

# Get feature data
featureInfo <- as.data.frame(mirExpr[,1:23])

save(featureInfo, file = "featureInfo.RData")

# Get sampleInfo
sampleInfo <- read.delim("sampleInfo.txt")
all(colnames(exprMatrix) == sampleInfo$SampleID)
save(sampleInfo, file = "sampleInfo.RData")





# Get normalized expression
exprMatrix <- as.matrix(mirExpr[,str_detect(colnames(mirExpr), "(norm)")])
exprMatrix <- apply(exprMatrix, 2, as.numeric)
rownames(exprMatrix) <- mirExpr$miRNA_Index
colnames(exprMatrix) <- str_remove_all(colnames(exprMatrix),"norm")
colnames(exprMatrix) <- str_remove_all(colnames(exprMatrix),"[()]")
exprMatrix <- log2(exprMatrix + 1)

