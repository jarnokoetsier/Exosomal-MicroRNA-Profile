#==============================================================================#

#   File:     PathwayAnalysis.R
#   Author:   Jarno Koetsier
#   Date:     February 6, 2023

#==============================================================================#


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
setwd("D:/RTTproject/ExosomeAnalysis/4. Pathway Analysis/")

# Get data directory
dataDir <- "D:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "D:/RTTproject/ExosomeAnalysis/4. Pathway Analysis/"


###############################################################################

# 1. Data collection and formatting

###############################################################################

# Load miRNA-Target Interaction (MTI) Network (miRtarbase v9.0)
MTI <- read_excel(paste0(dataDir,"hsa_MTI.xlsx"))

# Load feature information
load(paste0(dataDir,"featureInfo.RData"))

# Put miRNA name in correct format
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")

# Remove duplicated features
featureInfo <- featureInfo[!duplicated(featureInfo$CompleteName),]

# Filter MTI for measured miRNAs
MTI_filtered <- MTI[MTI$miRNA %in% featureInfo$CompleteName,]

# Save filtered MTI
save(MTI_filtered, file = paste0(outputDir,"MTI_filtered.RData"))

###############################################################################

# Data Preparation

###############################################################################

# Load mirRNA-target interaction data
load("MTI_filtered.RData")

#read GMT file
gmt <- clusterProfiler::read.gmt.wp("wikipathways-20220810-gmt-Homo_sapiens.gmt.txt")

# link pathway to gene
path2gene <- gmt[,c("wpid", "gene")]
path2name <- gmt[,c("wpid", "name")]


# Get the relevent gene sets (Run only the relevent secton):

#==============================================================================#
# Early expressed miRNAs
#==============================================================================#

# Get early expressed miRNAs
miRNA_set <- read.table("D:/RTTproject/ExosomeAnalysis/2. Time Analysis/DownExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]

gmt_fil <- gmt[gmt$wpid == "WP61",]
notch <- geneSet[geneSet$`Target Gene (Entrez ID)` %in% gmt_fil$gene,]
notch1 <- unique(notch[,c(2,4)])
table(notch1$`Target Gene`)

#==============================================================================#
# Lately expressed miRNAs
#==============================================================================#

# Get lately expressed miRNAs
miRNA_set <- read.table("D:/RTTproject/ExosomeAnalysis/2. Time Analysis/UpExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]

gmt_fil <- gmt[gmt$wpid == "WP4806",]
gmt_fil <- gmt[gmt$wpid == "WP2037",]
gmt_fil <- gmt[gmt$wpid == "WP3935",]
notch <- geneSet[geneSet$`Target Gene (Entrez ID)` %in% gmt_fil$gene,]
notch1 <- unique(notch[,c(2,4)])
table(notch1$`Target Gene`)

#==============================================================================#
# Ubiquitously expressed miRNAs
#==============================================================================#

# Get ubiquitously expressed miRNAs
miRNA_set <- read.table("D:/RTTproject/ExosomeAnalysis/2. Time Analysis/ConstantExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]

#==============================================================================#
# Differentially expressed miRNAs (RTT vs IC)
#==============================================================================#

# Get differentially expressed miRNAs (RTT vs IC
miRNA_set <- read.table("D:/RTTproject/ExosomeAnalysis/3. RTT vs IC/RTTvsIC.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNA_set$miRNA,]

gmt_fil <- gmt[gmt$wpid == "WP722",]
serotonin <- geneSet[geneSet$`Target Gene (Entrez ID)` %in% gmt_fil$gene,]
serotonin1 <- unique(serotonin[,c(2,4)])
table(serotonin1$`Target Gene`)

###############################################################################

# WikiPathays Overrepresentation Analysis

###############################################################################

# Perform ORA
WP <- enricher(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
               universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               TERM2GENE = path2gene,
               TERM2NAME = path2name)

result_WP <- WP@result
result_WP_copy <- WP@result[,c(1,2,5)]

# Number of permutations
nPerm <- 10000

# All miRNAs
all_mirnas <- unique(MTI_filtered$miRNA)

# Number of miRNAs to select during each permutation
n_miRNA <- length(miRNA_set$miRNA)

# Perform permutation analysis
set.seed(123)
for (i in 1:nPerm){
  
  # Get random miRNAs
  miRNA <- sample(all_mirnas, n_miRNA)
  
  # Get asssociated genes
  geneSet_perm <- MTI_filtered[MTI_filtered$miRNA %in% miRNA,]
  
  # Perform WP ORA
  WP <- enricher(gene = as.character(unique(geneSet_perm$`Target Gene (Entrez ID)`)),
                 universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                 pAdjustMethod = "fdr",
                 pvalueCutoff = Inf,
                 qvalueCutoff = Inf,
                 TERM2GENE = path2gene,
                 TERM2NAME = path2name)
  
  # Get results
  results <- WP@result[,c(1,5)]
  colnames(results) <- c("ID", paste0("Perm",i))
  
  # Combine results
  result_WP_copy <- left_join(result_WP_copy, results, by = c("ID" = "ID"))
}

# Get permutation results
perm_results <- -log10(result_WP_copy[,4:(nPerm+3)])
perm_results[is.na(perm_results)] <- 0

# Count number of times with more enrichment
test <- rowSums(perm_results > -log10(result_WP_copy[,3]))/nPerm

# Combine into data.frame
Output <- cbind.data.frame(result_WP_copy[,1:3], test)

# Calculate FDR
Output$FDR <- p.adjust(test, method = "fdr")

# Change column names
colnames(Output) <- c("ID", "Description", "pvalue (no perm)", "pvalue (perm)", "FDR")

# Order the rows by pvalue
Output <- arrange(Output, by = `pvalue (perm)`)

# Write output
save(perm_results, file = "perm_results_RTTvsIC.RData")
write.csv(Output, file = "WP_RTTvsIC.csv")
