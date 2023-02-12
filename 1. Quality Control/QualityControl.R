
#==============================================================================#

#   File:     QualityControl.R
#   Author:   Jarno Koetsier
#   Date:     February 6, 2023

#==============================================================================#


###############################################################################

# 0. Preparation

###############################################################################


# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# load packages
library(tidyverse)
library(edgeR)
library(heatmaply)
library(RColorBrewer)
library(ggpubr)
library(ggdendro)
library(ggsci)
library(grid)
library(ggpubr)
library(ggvenn)
library(GGally)
library(pcaMethods)
library(gridExtra)

# Set working directory
setwd("E:/RTTproject/ExosomeAnalysis/1. Quality Control")

# Get data directory
dataDir <- "E:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "E:/RTTproject/ExosomeAnalysis/1. Quality Control/"

###############################################################################

# 1. Data collection and formatting

###############################################################################

# Load data:

# 1. raw expression data
load(paste0(dataDir,"exprMatrix_raw.RData"))

# 2. feature information
load(paste0(dataDir,"featureInfo.RData"))

# 3. sample information
load(paste0(dataDir,"sampleInfo.RData"))

# Remove duplicated rows
exprMatrix <- exprMatrix[!duplicated(exprMatrix),]

# Put miRNA name in correct format
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")

# Remove duplicated features
featureInfo <- featureInfo[!duplicated(featureInfo$CompleteName),]


###############################################################################

# 2. Data Normalization

###############################################################################

# Collect samples
samples <- sampleInfo$SampleID

# Combine region/tissue (iPSC, Dorsal, Ventral), time (D0, D13, D40, and D75), and group (IC and RTT)
sampleInfo$Time_Group <- paste0(sampleInfo$Tissue, "_", sampleInfo$Time, "_", sampleInfo$Group)
time_group <- factor(sampleInfo$Time_Group)

# Make DGE list object
y <- DGEList(counts = exprMatrix[, samples],
             group = time_group)

# remove lowly expressed miRNAs
keep <- filterByExpr(y, min.count = 32)
y <- y[keep,,keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)

# Get normalized (log2) expression
logcpm <- log2(cpm(y)+1)

# Save normalized expression
save(logcpm, file = paste0(dataDir,"NormalizedExpr.RData"))


###############################################################################

# 3. Make QC plots

###############################################################################

# Format data for plotting
exprPlot <- gather(as.data.frame(logcpm))
exprPlot <- inner_join(exprPlot, sampleInfo, by = c("key" = "SampleID"))


#*****************************************************************************#
#   3.1. Boxplots
#*****************************************************************************#

# Make boxplot
exprBoxplot <- ggplot(exprPlot, aes(x = key, y = value, fill = Group)) +
  geom_boxplot(alpha = 0.8) +
  ylab(expression(log[2]~"normalized count")) +
  xlab("")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom", 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic")) +
  scale_fill_brewer(palette = "Set1")

# Save plot
ggsave(exprBoxplot, file = paste0(outputDir,"QC_Boxplots.png"), width = 8, height = 6)

#*****************************************************************************#
#   Density plots
#*****************************************************************************#

# Make density plot
exprDensityplot <- ggplot(exprPlot, aes(x = value, color = key)) +
  geom_density() +
  xlab(expression(log[2]~"normalized count")) +
  ylab("Density") +
  theme_classic() +
  theme_classic() +
  theme(axis.text.x = element_text(),
        legend.position = "none", 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic")) +
  scale_color_viridis_d()

# Save plot
ggsave(exprDensityplot, file = paste0(outputDir,"QC_Densityplots.png"), width = 8, height = 6)


#*****************************************************************************#
#   PCA plot
#*****************************************************************************#

# Perform PCA
pcaList <- pca(t(logcpm),
               scale = "none",
               center = TRUE,
               nPcs = 10)

# Calculate explained variance:
explVar <- round(pcaList@R2 * 100,2)

# Retrieve scores from pcaList object:
PCAscores <- as.data.frame(pcaList@scores)
PCAscores$SampleID <- rownames(PCAscores)
PCAscores <- inner_join(PCAscores, sampleInfo, by = c("SampleID" = "SampleID"))
PCAscores$Replicate <- as.character(PCAscores$Replicate)
PCAscores$Tissue[PCAscores$Tissue == "Cell"] <- "iPSC"

# Group: RTT vs IC
scorePlot_group <- ggplot(data = PCAscores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  ggtitle("Group") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_manual(breaks = c("RTT", "IC"),
                     values=c("#D61C4E", "#293462"))

# Time: D0, D13, D40, D75
scorePlot_time <- ggplot(data = PCAscores, aes(x = PC1, y = PC2, color = Time)) +
  geom_point(size = 2) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  ggtitle("Time") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_brewer(palette = "Set1")

# Region: Dorsal, Ventral, iPSC
scorePlot_tissue <- ggplot(data = PCAscores, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 2) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  ggtitle("Region") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_brewer(palette = "Dark2")

# Replicate
scorePlot_replicate <- ggplot(data = PCAscores, aes(x = PC1, y = PC2, color = Replicate)) +
  geom_point(size = 2) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  ggtitle("Replicate") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_viridis_d()

# Save plot
ggsave(paste0(outputDir, "PCAplot.png"), 
       grid.arrange(scorePlot_group, 
                    scorePlot_time, 
                    scorePlot_tissue,
                    scorePlot_replicate, 
                    ncol = 2, 
                    nrow = 2), 
       width = 8, height = 8)




