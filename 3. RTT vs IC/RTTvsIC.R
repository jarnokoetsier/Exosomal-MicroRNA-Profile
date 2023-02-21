#==============================================================================#

#   File:     RTTvsIC.R
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

# Set working directory
setwd("E:/RTTproject/ExosomeAnalysis/3. RTT vs IC/")

# Get data directory
dataDir <- "E:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "E:/RTTproject/ExosomeAnalysis/3. RTT vs IC/"


###############################################################################

# 1. Data collection and formatting

###############################################################################

# Load data:

# 1. raw expression data
load(paste0(dataDir,"exprMatrix_raw.RData"))

# 2. feature information
load(paste0(dataDir,"featureInfo.RData"))

# 3. sample information (IC only)
load(paste0(dataDir,"sampleInfo.RData"))

# Remove duplicated rows
exprMatrix <- exprMatrix[!duplicated(exprMatrix),]

# Put miRNA name in correct format
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")

# Remove duplicated features
featureInfo <- featureInfo[!duplicated(featureInfo$CompleteName),]

###############################################################################

# 2. Statistical analysis: RTT vs IC

###############################################################################

# Collect sample IDs
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

# Make design matrix
design <- model.matrix(~ 0 + time_group)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(y, design)

# Make Contrasts: RTT vs IC at every time point
all_contrasts <- makeContrasts(
  Cell_D0_RTTvsIC = time_groupCell_D0_RTT-time_groupCell_D0_IC,
  Dorsal_D13_RTTvsIC = time_groupDorsal_D13_RTT-time_groupDorsal_D13_IC,
  Dorsal_D40_RTTvsIC = time_groupDorsal_D40_RTT-time_groupDorsal_D40_IC,
  Dorsal_D75_RTTvsIC = time_groupDorsal_D75_RTT-time_groupDorsal_D75_IC,
  Ventral_D13_RTTvsIC = time_groupVentral_D13_RTT-time_groupVentral_D13_IC,
  Ventral_D40_RTTvsIC = time_groupVentral_D40_RTT-time_groupVentral_D40_IC,
  Ventral_D75_RTTvsIC = time_groupVentral_D75_RTT-time_groupVentral_D75_IC,
  levels = design
)

# Perform statistical analysis for each comparison and collect
# top table in a list object
topList <- list()
for (i in 1:7){
  # Test for significance
  test <- glmQLFTest(fit, contrast = all_contrasts[,i])
  
  #Get top table
  topTable <- topTags(test, n = Inf)$table
  topTable$ID <- rownames(topTable)
  
  # Collect top table into list
  topList[[i]] <- topTable
}

names(topList) <- colnames(all_contrasts)


# Collect the miRNAs that are significant for at least one time point
sigMir <- NULL
for (i in 1:length(topList)){
  sigMir <- c(sigMir,topList[[i]]$ID[topList[[i]]$FDR < 0.05])
}
sigMir <- unique(sigMir)
sigMir_name <- featureInfo$CompleteName[featureInfo$miRNA_Index %in% sigMir]

# Save significant miRNAs
write.table(sigMir_name, file = paste0(outputDir,"RTTvsIC.txt"), 
            quote = FALSE, col.names = FALSE, row.names = FALSE)


# Collect logFCs in a data frame
mir_logFC <- topList[[1]][,c(6,1)]
colnames(mir_logFC) <- c("ID", names(topList)[1])
for (i in 2:length(topList)){
  temp <- topList[[i]][,c(6,1)]
  colnames(temp) <- c("ID", names(topList)[i])
  mir_logFC <- inner_join(mir_logFC, temp, by = c("ID" = "ID"))
}


# Collext p-values in a data frame
mir_pvalue <- topList[[1]][,c(6,5)]
colnames(mir_pvalue) <- c("ID", names(topList)[1])
for (i in 2:length(topList)){
  temp <- topList[[i]][,c(6,5)]
  colnames(temp) <- c("ID", names(topList)[i])
  mir_pvalue <- inner_join(mir_pvalue, temp, by = c("ID" = "ID"))
}

###############################################################################

# 2. Prepare data for plotting

###############################################################################

# Put logFCs in correct format
plotData <- gather(mir_logFC[,2:8])
plotData$mir_index <- rep(mir_logFC$ID,7)
plotData$test <- paste0(plotData$key, plotData$mir_index)

# Put p-values in correct format
plotData_p <- gather(mir_pvalue[,2:8])
plotData_p$mir_index <- rep(mir_pvalue$ID,7)
plotData_p$test <- paste0(plotData_p$key, plotData_p$mir_index)
plotData_p <- plotData_p[,c(2,4)]
colnames(plotData_p) <- c("FDR", "test")

# Combine logFCs and p-values into a single data frame
plotData <- inner_join(plotData, plotData_p, by = c('test' = "test"))

# Add a column that indicates statistical signficance (FDR < 0.05)
plotData$Sig <- ifelse(plotData$FDR < 0.05, "Yes", "No")

# Filter for significant miRNAs only
plotData <- plotData[plotData$mir_index %in% sigMir,]

# Combine with feature info
plotData <- inner_join(plotData, featureInfo, by = c("mir_index" = "miRNA_Index"))

# Combine with sample info
plotData$key <- str_remove_all(plotData$key, "_RTTvsIC")
sampleInfo$tissue_time <- paste(sampleInfo$Tissue, sampleInfo$Time, sep = "_")
plotData <- inner_join(plotData, sampleInfo[(sampleInfo$Group == "IC") &
                                              (sampleInfo$Replicate == 1),], by = c("key" = "tissue_time"))

# Make sure that the groups are in the correct order for plotting
plotData$key <- factor(plotData$key,levels = c("Dorsal_D75",
                                              "Dorsal_D40",
                                              "Dorsal_D13",
                                              "Cell_D0",
                                              "Ventral_D13",
                                              "Ventral_D40",
                                               "Ventral_D75"))

# Change "Cell" to "iPSC
plotData$Tissue[plotData$Tissue == "Cell"] <- "iPSC"

# Make sure that the regions/tissues are in the correct order
plotData$Tissue <- factor(plotData$Tissue, levels = c("Dorsal", "iPSC", "Ventral"))

# The order of miRNAs is determined by clustering
rownames(mir_logFC) <- mir_logFC$ID
model <- hclust(dist(mir_logFC[,2:8]), "ward.D2")
rownames(featureInfo) <- featureInfo$miRNA_Index
mirna_ordered <- featureInfo[model$labels[model$order], "CompleteName"]
plotData$CompleteName <- factor(plotData$CompleteName, 
                                levels = mirna_ordered)

# Make main plot: Heatmap
mainPlot <- ggplot(data = plotData, aes(x = key, y = CompleteName, fill = value,  color = Sig)) +
  geom_tile(linewidth = 0.5, width = 0.9, height=0.7) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, trans = "pseudo_log") +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = "logFC") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())


# Row side color: Chr14 miRNA cluster (C14MC) vs others:

# Get miRNAs that belong to the C14MC:
cluster_chr14 <- featureInfo[featureInfo$genomeID == "chr14",]
cluster_chr14$start <- as.numeric(cluster_chr14$start)
cluster_chr14 <- cluster_chr14[(cluster_chr14$start < 101100000) &
                                 (cluster_chr14$start > 100800000),]
cluster_chr14 <- cluster_chr14[cluster_chr14$strand == "+",]

# Assign correct labels to the miRNAs
rowSideColor <- unique(data.frame(
  miRNA = plotData$CompleteName,
  Chr14 = ifelse(plotData$CompleteName %in% cluster_chr14$CompleteName, "C14MC", "Other")))

# Make row side colors
rowSideColorPlot <- ggplot() +
  geom_tile(data = rowSideColor, aes(x = miRNA, y = "label", fill = Chr14)) +
  coord_flip() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_manual(values = c("#F99417", "grey"))

# Column side color: time and tissue/region
colSideColor<- unique(data.frame(
  sample = plotData$key,
  time = plotData$Time,
  tissue = plotData$Tissue))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# Tissue/region
colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[(str_detect(colSideColor$sample, "D40") |
                                     str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots into a single figure
finalPlot <- ggarrange(NULL,
                       colSideColorPlot_tissue,
                       rowSideColorPlot,
                       mainPlot,
                       NULL,
                       colSideColorPlot_time,
                       heights = c(0.5,9,0.5), widths = c(2,8),nrow = 3,ncol = 2,
                       common.legend = FALSE)

# Save plot
ggsave(finalPlot, file = paste0(outputDir,"RTTvsIC_FinalPlot.png"), width = 8, height = 10)


# Get legends:

# Legend of logFCs in heatmap
legendPlot <- ggplot(data = plotData, aes(x = key, y = CompleteName, fill = value, color = Sig)) +
  geom_tile(linewidth = 0.5, width= 1, height=0.7) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       trans = "pseudo_log", breaks = c(-4,0,4)) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = expression(log[2]*FC)) +
  guides(color = "none") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank())


legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)


# Legend of row side color
legendPlot <- ggplot() +
  geom_tile(data = rowSideColor, aes(x = miRNA, y = "label", fill = Chr14)) +
  coord_flip() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_manual(values = c("#F99417", "grey"))
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)



###############################################################################

# 3. Plot Chromosome 14 cluster

###############################################################################


# Get members of the chr14 miRNA cluster (C14MC)
cluster_chr14 <- featureInfo[featureInfo$genomeID == "chr14",]
cluster_chr14$start <- as.numeric(cluster_chr14$start)
cluster_chr14 <- cluster_chr14[(cluster_chr14$start < 101100000) &
                                 (cluster_chr14$start > 100800000),]
cluster_chr14 <- cluster_chr14[cluster_chr14$strand == "+",]

# Get the expression of these C14MC miRNAs
load(paste0(dataDir,"NormalizedExpr.RData"))
expr_chr14 <- logcpm[rownames(logcpm) %in% cluster_chr14$miRNA_Index,]

# Put the expression matrix in correct format for plotting
df_chr14 <- gather(as.data.frame(expr_chr14), key = "SampleID", value = "Expr")
df_chr14$mir_id <- rep(rownames(expr_chr14), ncol(expr_chr14))

# Combine with feature information
df_chr14 <- inner_join(df_chr14, featureInfo, by = c("mir_id" = "miRNA_Index"))

# Combine with sample information
df_chr14 <- inner_join(df_chr14, sampleInfo, by = c("SampleID" = "SampleID"))

# Location of miRNA should be numerical
df_chr14$start <- as.numeric(df_chr14$start)

# Get sample groups
df_chr14$rep_group <- paste(df_chr14$Group, df_chr14$Tissue, df_chr14$Time, df_chr14$miR_name, sep = "_")

# Get average expression per group (IC, RTT), Time (D0-D75), and Region (Dorsal, ventral)
plotData <- df_chr14 %>%
  group_by(rep_group) %>%
  summarize(avgExpr = mean(Expr),
            start = start,
            Group = Group,
            Tissue = Tissue,
            Time = Time,
            SampleID = SampleID,
            miR_name = miR_name)

# Correct naming of sample groups (e.g., RTT: D40)
plotData$TimeGroup <- paste0(plotData$Group, ": ", plotData$Time)

# Plot ventral region
p_ventral <- ggplot() +
  geom_point(data = plotData[(plotData$Tissue != "Dorsal"),], 
             aes(x = start, y = avgExpr, color = TimeGroup, shape = Group)) +
  scale_color_manual(values = c(brewer.pal(n = 5, "Blues")[2:5],
                                brewer.pal(n = 5, "Reds")[2:5])) +
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  xlab("Location on chromosome 14") +
  ylab(expression("Mean "*log[2]* " expression")) +
  ylim(c(0,15)) +
  ggtitle("Ventral") +
  guides(shape = "none") +
  labs(shape = NULL, color = "Time") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        legend.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 12))
# Plot dorsal region
p_dorsal<- ggplot() +
  geom_point(data = plotData[(plotData$Tissue != "Ventral"),], 
             aes(x = start, y = avgExpr, color = TimeGroup, shape = Group)) +
  scale_color_manual(values = c(brewer.pal(n = 5, "Blues")[2:5],
                                brewer.pal(n = 5, "Reds")[2:5])) +
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  xlab(NULL) +
  ylab(expression("Mean "*log[2]* " expression")) +
  ylim(c(0,15)) +
  ggtitle("Dorsal") +
  guides(shape = "none") +
  labs(shape = NULL, color = "Time") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        legend.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 12)
  )

# Combine plots into single figure
p <- ggarrange(p_dorsal, p_ventral, nrow = 2, common.legend = TRUE, legend = "right")

# Save plot
ggsave(p, file = paste0(outputDir,"exprChr14Cluster.png"), height = 6, width = 8)

