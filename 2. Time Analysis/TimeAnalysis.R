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
setwd("E:/RTTproject/ExosomeAnalysis/2. Time Analysis")

# Get data directory
dataDir <- "E:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "E:/RTTproject/ExosomeAnalysis/2. Time Analysis"

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
save(mirList, file = file(outputDir,"mirList.RData"))

###############################################################################

# 5. Prepare data for plotting

###############################################################################

# Load significant miRNAs (if needed)
load("mirList.RData")

# Load normalized expression data
load(paste0(dataDir,"NormalizedExpr.RData"))
logcpm <- logcpm[,sampleInfo$SampleID]

# Adjust sample ID
sampleInfo$SampleID <- str_remove_all(sampleInfo$SampleID, "IC_MeCP2_R255X_")

# Label significant miRNAs
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

# Unit variance scale the expression matrix
exprMatrix_scaled <- (logcpm - rowMeans(logcpm))/(apply(logcpm,1,sd))
colnames(exprMatrix_scaled) <- str_remove_all(colnames(logcpm), "IC_MeCP2_R255X_")
rownames(exprMatrix_scaled) <- rownames(logcpm)

# Put expression matrix in correct format
exprDF <- gather(as.data.frame(exprMatrix_scaled), key = "SampleID", value = "Expr")
exprDF$miRNA_Index <- rep(rownames(exprMatrix_scaled), ncol(exprMatrix_scaled))

# Combine expression information with feature information
exprDF <- inner_join(exprDF, featureInfo[,c(2:9,24)], by = c("miRNA_Index" = "miRNA_Index"))

# Combine and filter for significant micrnas
exprDF_sig <- inner_join(exprDF, mirnaDF, by = c("CompleteName" = "mirName"))
length(unique(exprDF_sig$miR_name))

# Combine with sample Information
exprDF_sig <- inner_join(exprDF_sig, sampleInfo, by = c("SampleID" = "SampleID"))

# Adjust data frame for plotting
plotData <- exprDF_sig
plotData$SampleID <- factor(plotData$SampleID,
                              levels = c(rev(sampleInfo$SampleID[sampleInfo$Tissue == "Dorsal"]),
                                         sampleInfo$SampleID[sampleInfo$Tissue == "Cell"],
                                         sampleInfo$SampleID[sampleInfo$Tissue == "Ventral"]))

# save plotting data
#save(plotData, file = "plotData.RData")


###############################################################################

# Make plot

###############################################################################

# Load data
load("mirList.RData")

plotData$Tissue <- str_replace_all(plotData$Tissue, "Cell", "iPSC")

# Order the miRNAs based on Ward.D2 clustering
model <- hclust(dist(exprMatrix_scaled[unique(plotData$miRNA_Index),sampleInfo$SampleID]), "ward.D2")
rownames(featureInfo) <- featureInfo$miRNA_Index
mirna_ordered <- featureInfo[model$labels[model$order], "CompleteName"]

# Order the mirna's in plotData datafram
plotData$CompleteName <- factor(plotData$CompleteName, 
                                levels = mirna_ordered)
plotData$Tissue <- factor(plotData$Tissue, levels = c("Dorsal", "iPSC", "Ventral"))

# Save clusters
clusters <- cutree(model, k = 3)
save(clusters, file = "clusters.RData")

# Make dendrogram
dendroDF <- data.frame(Name = names(cutree(model, k = 3)), Cluster = cutree(model, k = 3))
#dendroDF$Cluster[dendroDF$Cluster == 1] <- 6L
#dendroDF$Cluster[dendroDF$Cluster == 2] <- 8L
#dendroDF$Cluster[dendroDF$Cluster == 3] <- 5L
#dendroDF$Cluster[dendroDF$Cluster == 4] <- 7L
dendroPlot <- ggplot() +
  geom_tile(data = dendroDF, 
            aes(x = 1, y = as.numeric(fct_reorder(Name, Cluster)), fill = factor(Cluster))) +
  geom_dendro(model,xlim = c(1.5,10), pointing = "side") +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = brewer.pal(n = 8, name = "Dark2")[c(4,6,8,7)])


# Make main plot
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



# Row side color
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

# Col side color
colSideColor<- unique(data.frame(
  sample = plotData$SampleID,
  time = plotData$Time,
  tissue = plotData$Tissue))


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

colSideColor_tissue <- unique(data.frame(
  sample = plotData$SampleID,
  tissue = plotData$Tissue))


colSideColorPlot_tissue <- ggplot(data = colSideColor_tissue) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[str_detect(colSideColor_tissue$sample, "2") & 
                                       (str_detect(colSideColor_tissue$sample, "D40") |
                                          str_detect(colSideColor_tissue$sample, "D0")),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots

# With dendrogram
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

ggsave(finalPlot, file = "TimeAnalysis_FinalPlot.png", width = 10, height = 8)

# Get legends
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

# constant expression

###############################################################################

# constant Expr
constantExpr <- logcpm[rowSums(logcpm > 12) == 21,]
mirna <- rownames(constantExpr)

write.table(featureInfo$CompleteName[featureInfo$miRNA_Index %in% mirna], file = "constantExpr.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)


# get expression
exprMatrix_filtered <- logcpm

# prepare data for plotting
df <- gather(as.data.frame(exprMatrix_filtered), key = "SampleID", value = "Expr")
df$mir_id <- rep(rownames(exprMatrix_filtered), ncol(exprMatrix_filtered))
df <- inner_join(df, featureInfo, by = c("mir_id" = "miRNA_Index"))
df<- inner_join(df, sampleInfo, by = c("SampleID" = "SampleID"))
df$start <- as.numeric(df$start)
df$rep_group <- paste(df$Group, df$Tissue, df$Time, df$miR_name, sep = "_")
df$Tissue <- str_replace_all(df$Tissue, "Cell", "iPSC")

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

plotData$plotGroup <- factor(paste(plotData$Tissue, plotData$Time,sep = "_"),
                             levels = c("Dorsal_D75",
                                        "Dorsal_D40",
                                        "Dorsal_D13",
                                        "iPSC_D0",
                                        "Ventral_D13",
                                        "Ventral_D40",
                                        "Ventral_D75"))


# Make plot
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

colSideColor<- unique(data.frame(
  sample = plotData$plotGroup,
  time = plotData$Time,
  tissue = plotData$Tissue))


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

p <- ggarrange(main, 
               colSideColorPlot_time,
               colSideColorPlot_tissue,
               nrow = 3,
               heights = c(8,1,1),
               align = "v",
               common.legend = TRUE,
               legend = "right")

# Save plot
ggsave(p, file = "constExpr_color_parallel.png", width = 10, height = 6)

###############################################################################

# 302-family

###############################################################################

# get expression
exprMatrix_filtered <- logcpm

# prepare data for plotting
df <- gather(as.data.frame(exprMatrix_filtered), key = "SampleID", value = "Expr")
df$mir_id <- rep(rownames(exprMatrix_filtered), ncol(exprMatrix_filtered))
df <- inner_join(df, featureInfo, by = c("mir_id" = "miRNA_Index"))
df<- inner_join(df, sampleInfo, by = c("SampleID" = "SampleID"))
df$start <- as.numeric(df$start)
df$rep_group <- paste(df$Group, df$Tissue, df$Time, df$miR_name, sep = "_")
df$Tissue <- str_replace_all(df$Tissue, "Cell", "iPSC")

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

plotData$plotGroup <- factor(paste(plotData$Tissue, plotData$Time,sep = "_"),
                             levels = c("Dorsal_D75",
                                        "Dorsal_D40",
                                        "Dorsal_D13",
                                        "iPSC_D0",
                                        "Ventral_D13",
                                        "Ventral_D40",
                                        "Ventral_D75"))


# Get 302-family
mirna <- featureInfo$miRNA_Index[(str_detect(featureInfo$CompleteName, "302")) | (str_detect(featureInfo$CompleteName, "367-"))]

# Make plot
main <- ggplot() +
  geom_line(data = plotData, aes(x = plotGroup, y = avgExpr, group = miR_name), color = "grey") +
  geom_line(data = plotData[plotData$miR_Index %in% mirna,], 
            aes(x = plotGroup, y = avgExpr, group = miR_name, color = miR_name), size = 1) +
  ylab(expression("Mean" ~ log[2] ~ "expression")) +
  labs(color = NULL) +
  scale_color_manual(values = c(brewer.pal(n = 7, name = "Dark2"), "red", "black")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

colSideColor<- unique(data.frame(
  sample = plotData$plotGroup,
  time = plotData$Time,
  tissue = plotData$Tissue))


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

p <- ggarrange(main, 
          colSideColorPlot_time,
          colSideColorPlot_tissue,
          nrow = 3,
          heights = c(8,1,1),
          align = "v",
          common.legend = TRUE,
          legend = "right")

# Save plot
ggsave(p, file = "mir302_color_parallel.png", width = 10, height = 6)


###############################################################################

# Up/Down over time

###############################################################################

# Sig miRNAs over time
load("clusters.RData")
mirna_yellow <- names(clusters[clusters == 2])
write.table(featureInfo$CompleteName[featureInfo$miRNA_Index %in% mirna_yellow], file = "DownExpr.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)


main <- ggplot() +
  geom_line(data = plotData, aes(x = plotGroup, y = avgExpr, group = miR_name), color = "grey") +
  geom_line(data = plotData[plotData$miR_Index %in% mirna_yellow,], 
            aes(x = plotGroup, y = avgExpr, group = miR_name), 
            color = brewer.pal(n = 8, name = "Dark2")[c(6)], size = 1) +
  ylab(expression("Mean "~log[2]~" expression")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

colSideColor<- unique(data.frame(
  sample = plotData$plotGroup,
  time = plotData$Time,
  tissue = plotData$Tissue))


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

p <- ggarrange(main, 
               colSideColorPlot_time,
               colSideColorPlot_tissue,
               nrow = 3,
               heights = c(8,1,1),
               align = "v")

ggsave(p, file = "DownOverTime_parallel.png", width = 8, height = 6)

mirna_pink <- names(clusters[clusters == 1])
write.table(featureInfo$CompleteName[featureInfo$miRNA_Index %in% mirna_pink], file = "UpExpr.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

main <- ggplot() +
  geom_line(data = plotData, aes(x = plotGroup, y = avgExpr, group = miR_name), color = "grey") +
  geom_line(data = plotData[plotData$miR_Index %in% mirna_pink,], 
            aes(x = plotGroup, y = avgExpr, group = miR_name), 
            color = brewer.pal(n = 8, name = "Dark2")[c(4)], size = 1) +
  ylab(expression("Mean "~log[2]~" expression")) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

colSideColor<- unique(data.frame(
  sample = plotData$plotGroup,
  time = plotData$Time,
  tissue = plotData$Tissue))


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

p <- ggarrange(main, 
               colSideColorPlot_time,
               colSideColorPlot_tissue,
               nrow = 3,
               heights = c(8,1,1),
               align = "v")

ggsave(p, file = "UpOverTime_parallel.png", width = 8, height = 6)


