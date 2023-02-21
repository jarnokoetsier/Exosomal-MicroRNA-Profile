#==============================================================================#

#   File:     TimeAnalysis.R
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
library(ggdendroplot)
library(ggsci)
library(grid)
library(ggpubr)
library(ggvenn)
library(GGally)

# Set working directory
setwd("E:/RTTproject/ExosomeAnalysis/2. Time Analysis/")

# Get data directory
dataDir <- "E:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "E:/RTTproject/ExosomeAnalysis/2. Time Analysis/"

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
sampleInfo <- sampleInfo[sampleInfo$Group == "IC",]
exprMatrix <- exprMatrix[,sampleInfo$SampleID]

# Remove duplicated rows
exprMatrix <- exprMatrix[!duplicated(exprMatrix),]

# Put miRNA name in correct format
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")

# Remove duplicated features
featureInfo <- featureInfo[!duplicated(featureInfo$CompleteName),]

# Make empty list (will be used to save significant miRNAs)
mirList <- list()

###############################################################################

# 2. Data Normalization and Fitting of Linear Model

###############################################################################

# Get samples from dorsal and ventral region
samples_ventral <- sampleInfo$SampleID[(sampleInfo$Tissue != "Dorsal")]
samples_dorsal <- sampleInfo$SampleID[(sampleInfo$Tissue != "Ventral")]

# Combine region/tissue (iPSC, Dorsal, Ventral), time (D0, D13, D40, and D75), and group (IC and RTT)
sampleInfo$Time_Group <- paste0(sampleInfo$Tissue, "_", sampleInfo$Time, "_", sampleInfo$Group)
time_group <- factor(sampleInfo$Time_Group)

# Make DGE list object
y <- DGEList(counts = exprMatrix,
             group = time_group)

# remove lowly expressed miRNAs
keep <- filterByExpr(y, min.count = 32)
y <- y[keep,,keep.lib.sizes=FALSE]

# Normalize the data
y <- calcNormFactors(y)

# Make design matrix
design <- model.matrix(~0 + time_group)

# Estimate dispersion
y <- estimateDisp(y, design)

# Fit linear model (quasi-likelihood F-tests)
fit <- glmQLFit(y, design)

###############################################################################

# 3. Time Comparison for Ventral Region (IC)

###############################################################################

# Make contrasts (every time point comparison possible)
all_contrasts <- makeContrasts(
  D13vsD0 = time_groupVentral_D13_IC - time_groupCell_D0_IC,
  D40vsD0 = time_groupVentral_D40_IC - time_groupCell_D0_IC,
  D75vsD0  = time_groupVentral_D75_IC - time_groupCell_D0_IC,
  D40vsD13 = time_groupVentral_D40_IC - time_groupVentral_D13_IC,
  D75vsD13 = time_groupVentral_D75_IC - time_groupVentral_D13_IC,
  D75vsD40 = time_groupVentral_D75_IC - time_groupVentral_D40_IC,
  levels = design
)

# Test for significance
test <- glmQLFTest(fit, contrast = all_contrasts)

# Make top table
topTable <- topTags(test, n = Inf)$table
topTable$ID <- rownames(topTable)

# Combine with feature info
topTable_combined <- inner_join(topTable, featureInfo, by = c("ID" = "miRNA_Index"))
rownames(topTable_combined) <- topTable_combined$ID

# Get significant miRNAs (FDR < 0.05)
sigMir <- topTable$ID[topTable$FDR < 0.05]

# Save significant miRNAs
mirList[["Ventral_IC"]] <- topTable_combined[sigMir, "CompleteName"]


###############################################################################

# 4. Time Comparison for Dorsal Region (IC)

###############################################################################

# Make contrasts (every time point comparison possible)
all_contrasts <- makeContrasts(
  D13vsD0 = time_groupDorsal_D13_IC - time_groupCell_D0_IC,
  D40vsD0 = time_groupDorsal_D40_IC - time_groupCell_D0_IC,
  D75vsD0  = time_groupDorsal_D75_IC - time_groupCell_D0_IC,
  D40vsD13 = time_groupDorsal_D40_IC - time_groupDorsal_D13_IC,
  D75vsD13 = time_groupDorsal_D75_IC - time_groupDorsal_D13_IC,
  D75vsD40 = time_groupDorsal_D75_IC - time_groupDorsal_D40_IC,
  levels = design
)

# Test for significance
test <- glmQLFTest(fit, contrast = all_contrasts)

# Make top table
topTable <- topTags(test, n = Inf)$table
topTable$ID <- rownames(topTable)

# Combine with feature info
topTable_combined <- inner_join(topTable, featureInfo, by = c("ID" = "miRNA_Index"))
rownames(topTable_combined) <- topTable_combined$ID

# Get significant miRNAs
sigMir <- topTable$ID[topTable$FDR < 0.05]

# Save significant miRNAs
mirList[["Dorsal_IC"]] <- topTable_combined[sigMir, "CompleteName"]

# Save mirList object
save(mirList, file = paste0(outputDir,"mirList.RData"))

###############################################################################

# 5. Prepare data for plotting

###############################################################################

# Load significant miRNAs (if needed)
load(paste0(outputDir,"mirList.RData"))

# Load normalized expression data
load(paste0(dataDir,"NormalizedExpr.RData"))

# Get IC samples only
logcpm <- logcpm[,sampleInfo$SampleID]

# Adjust sample ID
sampleInfo$SampleID <- str_remove_all(sampleInfo$SampleID, "IC_MeCP2_R255X_")

# Significant miRNAs: Ventral only, Dorsal only, or both
both_sig <- intersect(mirList$Ventral_IC, mirList$Dorsal_IC)
ventral_sig <- setdiff(mirList$Ventral_IC, mirList$Dorsal_IC)
dorsal_sig <- setdiff(mirList$Dorsal_IC, mirList$Ventral_IC)

# Put significance labels of miRNAs in data frame
mirnaDF <- data.frame(mirName = c(both_sig,
                                  ventral_sig, 
                                  dorsal_sig),
                      Sig = c(rep("Both", length(both_sig)),
                              rep("Ventral", length(ventral_sig)),
                              rep("Dorsal", length(dorsal_sig)))
                      )

# Unit variance scale the expression matrix: (x - mean(x))/sd(x)
exprMatrix_scaled <- (logcpm - rowMeans(logcpm))/(apply(logcpm,1,sd))
colnames(exprMatrix_scaled) <- str_remove_all(colnames(logcpm), "IC_MeCP2_R255X_")
rownames(exprMatrix_scaled) <- rownames(logcpm)

# Put expression matrix in correct format for plotting with ggplot
exprDF <- gather(as.data.frame(exprMatrix_scaled), key = "SampleID", value = "Expr")
exprDF$miRNA_Index <- rep(rownames(exprMatrix_scaled), ncol(exprMatrix_scaled))

# Combine expression information with feature information
exprDF <- inner_join(exprDF, featureInfo[,c(2:9,24)], by = c("miRNA_Index" = "miRNA_Index"))

# Combine and filter for significant miRNAs
exprDF_sig <- inner_join(exprDF, mirnaDF, by = c("CompleteName" = "mirName"))
length(unique(exprDF_sig$miR_name))

# Combine with sample Information
exprDF_sig <- inner_join(exprDF_sig, sampleInfo, by = c("SampleID" = "SampleID"))

# Convert sample IDs to a factor to ensure correct order of samples in plot
plotData <- exprDF_sig
plotData$SampleID <- factor(plotData$SampleID,
                              levels = c(rev(sampleInfo$SampleID[sampleInfo$Tissue == "Dorsal"]),
                                         sampleInfo$SampleID[sampleInfo$Tissue == "Cell"],
                                         sampleInfo$SampleID[sampleInfo$Tissue == "Ventral"]))


# Change "Cell" to "iPSC"
plotData$Tissue <- str_replace_all(plotData$Tissue, "Cell", "iPSC")


###############################################################################

# 6. Cluster the miRNAs and make dendrogram

###############################################################################

# Perform hierarchical clustering
model <- hclust(dist(exprMatrix_scaled[unique(plotData$miRNA_Index),sampleInfo$SampleID],
                     method = "euclidean"), "ward.D2")

# Make clusters
clusters <- cutree(model, k = 3)

# Get miRNA order (needed for plotting)
rownames(featureInfo) <- featureInfo$miRNA_Index
mirna_ordered <- featureInfo[model$labels[model$order], "CompleteName"]

# Put cluster information into data frame for plotting
dendroDF <- data.frame(Name = factor(mirna_ordered,
                                     levels = mirna_ordered), 
                       Cluster = cutree(model, k = 3))

# Make dendrogram
dendroPlot <- ggplot() +
  geom_tile(data = dendroDF, 
            aes(x = 1, y = as.numeric(fct_reorder(Name, Cluster)), fill = factor(Cluster))) +
  geom_dendro(model,xlim = c(1.5,10), pointing = "side") +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = brewer.pal(n = 8, name = "Dark2")[c(4,6,8)])



###############################################################################

# 7. Make final plot

###############################################################################

# Order the miRNAs in the data frame
plotData$CompleteName <- factor(plotData$CompleteName, 
                                levels = mirna_ordered)
plotData$Tissue <- factor(plotData$Tissue, levels = c("Dorsal", "iPSC", "Ventral"))

# Make heatmap of expression values (= main plot)
mainPlot <- ggplot(data = plotData, aes(x = SampleID, y = CompleteName, fill = Expr)) +
  geom_tile(color = "white") +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  scale_fill_viridis_c(limits = c(-2, 2), oob = scales::squish) +
  labs(fill = expression("Standardized\n" ~log[2] ~ "expression")) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())


# Add row side color: significant in dorsal, ventral, or both regions?
rowSideColor <- unique(data.frame(
  miRNA = plotData$CompleteName,
  Sig = factor(plotData$Sig, levels = c("Dorsal", "Ventral", "Both"))))

rowSideColorPlot <- ggplot() +
  geom_tile(data = rowSideColor, aes(x = miRNA, y = "label", fill = Sig)) +
  coord_flip() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

# Add row side color: region/tissue and time
colSideColor<- unique(data.frame(
  sample = plotData$SampleID,
  time = plotData$Time,
  tissue = plotData$Tissue))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor[str_detect(colSideColor$sample, "2"),],
            aes(x = sample, y = "label", label = time)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),#element_text(angle = 90),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Reds")

# Region/Tissue
colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[str_detect(colSideColor$sample, "2") & 
                                       (str_detect(colSideColor$sample, "D40") |
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
                       NULL,
                       rowSideColorPlot,
                       mainPlot,
                       dendroPlot,
                       NULL,
                       colSideColorPlot_time,
                       NULL,
                       heights = c(0.5,9,0.5), widths = c(2,8,2),nrow = 3,ncol = 3,
                       common.legend = FALSE)

# Save plot
ggsave(finalPlot, file = paste0(outputDir,"TimeAnalysis_FinalPlot.png"), width = 10, height = 8)


# Get legends for plots:

# Legend of expression values
legendPlot <- ggplot() +
  geom_tile(data = plotData, aes(x = SampleID, y = CompleteName, fill = Expr)) +
  scale_fill_viridis_c(limits = c(-2, 2), oob = scales::squish) +
  labs(fill = expression("Standardized\n" ~log[2] ~ "expression")) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right")
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)

# Legend of row labels
legendPlot <-ggplot() +
  geom_tile(data = rowSideColor, aes(x = miRNA, y = "label", fill = Sig)) +
  labs(fill = "FDR-adjusted\np-value < 0.05") +
  coord_flip() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_brewer(palette = "Set1")
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)


###############################################################################

# 8. Identification of ubiquitously expressed miRNAs

###############################################################################

# miRNAs with a constantly high expression (ubiquitously expressed):
# expression higher than 12 in all 21 isogenic control samples
constantExpr <- logcpm[rowSums(logcpm > 12) == 21,]
mirna <- rownames(constantExpr)

# Save the miRNAs in a file (will be used for pathway analysis)
write.table(featureInfo$CompleteName[featureInfo$miRNA_Index %in% mirna], file = paste0(outputDir,"constantExpr.txt"), 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Prepare data for plotting
exprMatrix_filtered <- logcpm
df <- gather(as.data.frame(exprMatrix_filtered), key = "SampleID", value = "Expr")
df$mir_id <- rep(rownames(exprMatrix_filtered), ncol(exprMatrix_filtered))

# Combine with feature information
df <- inner_join(df, featureInfo, by = c("mir_id" = "miRNA_Index"))

# Combine with sample information
load(paste0(dataDir,"sampleInfo.RData"))
df<- inner_join(df, sampleInfo, by = c("SampleID" = "SampleID"))

# Format data
df$start <- as.numeric(df$start)
df$rep_group <- paste(df$Group, df$Tissue, df$Time, df$miR_name, sep = "_")
df$Tissue <- str_replace_all(df$Tissue, "Cell", "iPSC")

# Get mean expression per time and region/tissue
plotData <- df %>%
  group_by(rep_group) %>%
  summarize(avgExpr = mean(Expr),
            start = start,
            Group = Group,
            Tissue = Tissue,
            Time = Time,
            SampleID = SampleID,
            miR_name = CompleteName,
            miR_Index = mir_id)

# Set order of groups for plotting
plotData$plotGroup <- factor(paste(plotData$Tissue, plotData$Time,sep = "_"),
                             levels = c("Dorsal_D75",
                                        "Dorsal_D40",
                                        "Dorsal_D13",
                                        "iPSC_D0",
                                        "Ventral_D13",
                                        "Ventral_D40",
                                        "Ventral_D75"))


# Make main plot
main <- ggplot() +
  geom_line(data = plotData, aes(x = plotGroup, y = avgExpr, group = miR_name), color = "grey") +
  geom_line(data = plotData[plotData$miR_Index %in% mirna,], 
            aes(x = plotGroup, y = avgExpr, group = miR_name, color = miR_name), size = 1) +
  ylab(expression("Mean" ~ log[2] ~ "expression")) +
  labs(color = NULL) +
  scale_color_manual(values = c(brewer.pal(n = 10, name = "Paired"))) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())


# Add column side colors
colSideColor<- unique(data.frame(
  sample = plotData$plotGroup,
  time = plotData$Time,
  tissue = plotData$Tissue))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
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
  geom_text(data = colSideColor[(colSideColor$time == "D40") |
                                  (colSideColor$time == "D0"),],
            aes(x = sample, y = "label", label = tissue)) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")

# Combine plots into single image
p <- ggarrange(main, 
               colSideColorPlot_time,
               colSideColorPlot_tissue,
               nrow = 3,
               heights = c(8,1,1),
               align = "v",
               common.legend = TRUE,
               legend = "right")

# Save plot
ggsave(p, file = paste0(outputDir, "constExpr_color_parallel.png"), width = 10, height = 6)



###############################################################################

# 9. Plot expression of hsa-miR-302/367 cluster

###############################################################################

# Prepare data for plotting
exprMatrix_filtered <- logcpm
df <- gather(as.data.frame(exprMatrix_filtered), key = "SampleID", value = "Expr")
df$mir_id <- rep(rownames(exprMatrix_filtered), ncol(exprMatrix_filtered))

# Combine with feature information
df <- inner_join(df, featureInfo, by = c("mir_id" = "miRNA_Index"))

# Combine with sample information
load(paste0(dataDir,"sampleInfo.RData"))
df<- inner_join(df, sampleInfo, by = c("SampleID" = "SampleID"))

# Format data
df$start <- as.numeric(df$start)
df$rep_group <- paste(df$Group, df$Tissue, df$Time, df$miR_name, sep = "_")
df$Tissue <- str_replace_all(df$Tissue, "Cell", "iPSC")

# Get mean expression per time and region/tissue
plotData <- df %>%
  group_by(rep_group) %>%
  summarize(avgExpr = mean(Expr),
            start = start,
            Group = Group,
            Tissue = Tissue,
            Time = Time,
            SampleID = SampleID,
            miR_name = CompleteName,
            miR_Index = mir_id)

# Set order of groups for plotting
plotData$plotGroup <- factor(paste(plotData$Tissue, plotData$Time,sep = "_"),
                             levels = c("Dorsal_D75",
                                        "Dorsal_D40",
                                        "Dorsal_D13",
                                        "iPSC_D0",
                                        "Ventral_D13",
                                        "Ventral_D40",
                                        "Ventral_D75"))


# Get hsa-miR-302/367 cluster
mirna <- featureInfo$miRNA_Index[(str_detect(featureInfo$CompleteName, "302")) | (str_detect(featureInfo$CompleteName, "367-"))]

# Main plot
main <- ggplot() +
  geom_line(data = plotData, aes(x = plotGroup, y = avgExpr, group = miR_name), color = "grey") +
  geom_line(data = plotData[plotData$miR_Index %in% mirna,], 
            aes(x = plotGroup, y = avgExpr, group = miR_name, color = miR_name), linewidth = 1) +
  ylab(expression("Mean" ~ log[2] ~ "expression")) +
  labs(color = NULL) +
  scale_color_manual(values = c(brewer.pal(n = 7, name = "Dark2"), "red", "black")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

# Column side colors
colSideColor<- unique(data.frame(
  sample = plotData$plotGroup,
  time = plotData$Time,
  tissue = plotData$Tissue))

# Time
colSideColorPlot_time <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = time)) +
  geom_text(data = colSideColor,
            aes(x = sample, y = "label", label = time)) +
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
  geom_text(data = colSideColor[(colSideColor$time == "D40") |
                                  (colSideColor$time == "D0"),],
            aes(x = sample, y = "label", label = tissue)) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")

# Combine plots into single figure
p <- ggarrange(main, 
          colSideColorPlot_time,
          colSideColorPlot_tissue,
          nrow = 3,
          heights = c(8,1,1),
          align = "v",
          common.legend = TRUE,
          legend = "right")

# Save plot
ggsave(p, file = paste0(outputDir,"mir302_color_parallel.png"), width = 10, height = 6)

