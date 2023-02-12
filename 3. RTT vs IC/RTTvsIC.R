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
setwd("E:/RTTproject/ExosomeAnalysis/RTTvsIC")

# Set homeDir
homeDir <- "E:/RTTproject/ExosomeAnalysis/"

###############################################################################

# Data collection

###############################################################################

# Load data
load(paste0(homeDir,"exprMatrix_raw.RData"))
load(paste0(homeDir,"featureInfo.RData"))
load(paste0(homeDir,"sampleInfo.RData"))

# remove duplicated rows
# MiRNAs with same name (but different chromosomes) have the same expression
# in the expression matrix
exprMatrix <- exprMatrix[!duplicated(exprMatrix),]

# Make name of samples
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")
featureInfo <- featureInfo[!duplicated(featureInfo$CompleteName),]

###############################################################################

# 1. IC vs RTT

###############################################################################

# Make DGE list object
samples <- sampleInfo$SampleID
# Combine region/tissue (iPSC, Dorsal, Ventral), time (D0, D13, D40, and D75), and group (IC and RTT)
sampleInfo$Time_Group <- paste0(sampleInfo$Tissue, "_", sampleInfo$Time, "_", sampleInfo$Group)
time_group <- factor(sampleInfo$Time_Group)

y <- DGEList(counts = exprMatrix[, samples],
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

# Contrasts
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

# Perform statistical analysis
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

# miRNA's that are significant at at least one time point
sigMir <- NULL
for (i in 1:length(topList)){
  sigMir <- c(sigMir,topList[[i]]$ID[topList[[i]]$FDR < 0.05])
}
sigMir <- unique(sigMir)
sigMir_name <- featureInfo$CompleteName[featureInfo$miRNA_Index %in% sigMir]

# save sig miRNAs
write.table(sigMir_name, file = "RTTvsIC.txt", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# Get logFC dataframe
mir_logFC <- topList[[1]][,c(6,1)]
colnames(mir_logFC) <- c("ID", names(topList)[1])
for (i in 2:length(topList)){
  temp <- topList[[i]][,c(6,1)]
  colnames(temp) <- c("ID", names(topList)[i])
  mir_logFC <- inner_join(mir_logFC, temp, by = c("ID" = "ID"))
}


# Get p-value dataframe
mir_pvalue <- topList[[1]][,c(6,5)]
colnames(mir_pvalue) <- c("ID", names(topList)[1])
for (i in 2:length(topList)){
  temp <- topList[[i]][,c(6,5)]
  colnames(temp) <- c("ID", names(topList)[i])
  mir_pvalue <- inner_join(mir_pvalue, temp, by = c("ID" = "ID"))
}

###############################################################################

# 2. Make plot

###############################################################################

# Collect logFCs
plotData <- gather(mir_logFC[,2:8])
plotData$mir_index <- rep(mir_logFC$ID,7)
plotData$test <- paste0(plotData$key, plotData$mir_index)

# Collect p-values
plotData_p <- gather(mir_pvalue[,2:8])
plotData_p$mir_index <- rep(mir_pvalue$ID,7)
plotData_p$test <- paste0(plotData_p$key, plotData_p$mir_index)
plotData_p <- plotData_p[,c(2,4)]
colnames(plotData_p) <- c("FDR", "test")

# Combine logFCs and pvalues
plotData <- inner_join(plotData, plotData_p, by = c('test' = "test"))
plotData$Sig <- ifelse(plotData$FDR < 0.05, "Yes", "No")

# Filter for significant miRNAs only
plotData <- plotData[plotData$mir_index %in% sigMir,]


# combine with feature info
plotData <- inner_join(plotData, featureInfo, by = c("mir_index" = "miRNA_Index"))

# combine with sample info
plotData$key <- str_remove_all(plotData$key, "_RTTvsIC")

# prepare data for plotting
sampleInfo$tissue_time <- paste(sampleInfo$Tissue, sampleInfo$Time, sep = "_")
plotData <- inner_join(plotData, sampleInfo[(sampleInfo$Group == "IC") &
                                              (sampleInfo$Replicate == 1),], by = c("key" = "tissue_time"))

plotData$key <- factor(plotData$key,levels = c("Dorsal_D75",
                                              "Dorsal_D40",
                                              "Dorsal_D13",
                                              "Cell_D0",
                                              "Ventral_D13",
                                              "Ventral_D40",
                                               "Ventral_D75"))

plotData$Tissue[plotData$Tissue == "Cell"] <- "iPSC"
plotData$Tissue <- factor(plotData$Tissue, levels = c("Dorsal", "iPSC", "Ventral"))


# Order of miRNAs
rownames(mir_logFC) <- mir_logFC$ID
model <- hclust(dist(mir_logFC[,2:8]), "ward.D2")
rownames(featureInfo) <- featureInfo$miRNA_Index
mirna_ordered <- featureInfo[model$labels[model$order], "CompleteName"]

# Order the mirna's in plotData datafram
plotData$CompleteName <- factor(plotData$CompleteName, 
                                levels = mirna_ordered)

# Main plot: Heatmap
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

# Row side color: Chr14 cluster vs others
cluster_chr14 <- featureInfo[featureInfo$genomeID == "chr14",]
cluster_chr14$start <- as.numeric(cluster_chr14$start)
cluster_chr14 <- cluster_chr14[(cluster_chr14$start < 101100000) &
                                 (cluster_chr14$start > 100800000),]
rowSideColor <- unique(data.frame(
  miRNA = plotData$CompleteName,
  Chr14 = ifelse(plotData$CompleteName %in% cluster_chr14$CompleteName, "C14MC", "Other")))


rowSideColorPlot <- ggplot() +
  geom_tile(data = rowSideColor, aes(x = miRNA, y = "label", fill = Chr14)) +
  coord_flip() +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(hjust = 1)) +
  #scale_fill_brewer(palette = "Set1") +
  #scale_fill_manual(values = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[c(4,6)])
  scale_fill_manual(values = c("#F99417", "grey"))

# Col side color: time and tissue
colSideColor<- unique(data.frame(
  sample = plotData$key,
  time = plotData$Time,
  tissue = plotData$Tissue))


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


# Combine plots
finalPlot <- ggarrange(NULL,
                       colSideColorPlot_tissue,
                       rowSideColorPlot,
                       mainPlot,
                       NULL,
                       colSideColorPlot_time,
                       heights = c(0.5,9,0.5), widths = c(2,8),nrow = 3,ncol = 2,
                       common.legend = FALSE)

# Save plot
ggsave(finalPlot, file = "RTTvsIC_FinalPlot.png", width = 8, height = 10)




# Get legends
legendPlot <- ggplot(data = plotData, aes(x = key, y = CompleteName, fill = value, color = Sig)) +
  geom_tile(linewidth = 0.5, width= 1, height=0.7) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0, 
                       trans = "pseudo_log", breaks = c(-4,0,4)) +
  scale_color_manual(values = c("white", "black")) +
  labs(fill = expression(log[2]*FC)) +
  guides(color = FALSE) +
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

# Parallel plot

###############################################################################

# Prepare data (logFCs) for plotting
plotData <- gather(mir_logFC[,2:8])
plotData$ID <- rep(mir_logFC$ID, 7)
plotData$Tissue <- str_remove_all(plotData$key, "_.*")
plotData$Time <- str_remove_all(plotData$key, "_RTTvsIC")
plotData$Time <- str_remove_all(plotData$Time, ".*_")


plotData$Tissue[plotData$Time == "D0"] <- "iPSC"
plotData$Sample <- factor(paste(plotData$Tissue, plotData$Time, sep = "_"),
                          levels = c("Dorsal_D75",
                                     "Dorsal_D40",
                                     "Dorsal_D13",
                                     "iPSC_D0",
                                     "Ventral_D13",
                                     "Ventral_D40",
                                     "Ventral_D75"))

plotData <- inner_join(plotData, featureInfo, by = c("ID" = "miRNA_Index"))


# Get miRNAs that are part of chr14 cluster
cluster_chr14 <- featureInfo[featureInfo$genomeID == "chr14",]
cluster_chr14$start <- as.numeric(cluster_chr14$start)
cluster_chr14 <- cluster_chr14[(cluster_chr14$start < 101100000) &
                                 (cluster_chr14$start > 100800000),]
cluster_chr14 <- cluster_chr14[cluster_chr14$strand == "+",]
mirna <- cluster_chr14$miRNA_Index


# Main plot: logFC over time
main <- ggplot() +
  geom_line(data = plotData, aes(x = Sample, y = value, group = ID), color = "grey") +
  geom_line(data = plotData[plotData$ID %in% mirna,], 
            aes(x = Sample, y = value, group = ID, color = reorder(start, as.numeric(start))), size = 1) +
  ylab("logFC") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_color_viridis_d(option = "magma")

# Side colors: tissue and time
colSideColor<- unique(data.frame(
  sample = plotData$Sample,
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

# Combine plots into single image
p <- ggarrange(main, 
               colSideColorPlot_time,
               colSideColorPlot_tissue,
               nrow = 3,
               heights = c(8,1,1),
               align = "v")

# Save output
ggsave(p, file = "chr14_parallel.png", width = 8, height = 6)


# Make legend
test <- plotData[!duplicated(plotData$start),]
p <- ggplot(test[test$ID %in% mirna,]) +
  geom_tile(aes(x = as.numeric(start), y = 1, fill = reorder(start, as.numeric(start)))) +
  theme_classic() +
  xlab("Location on chromosome 14") +
  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  coord_flip() +
  scale_fill_viridis_d(option = "magma") +
  theme(legend.position = "none",
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 1),
        axis.text.y = element_text(angle = 90, hjust = 0),
        axis.ticks.y = element_line(size = 2))

# Save legend
ggsave(p, file = "legendParallel.png", width = 1, height = 8)

























plotData_ventral <- inner_join(topTable_ventral, featureInfo, by = c("ID" = "miRNA_Index"))
colnames(plotData_ventral) <- c("D0", "D13", "D40", "D75", colnames(plotData_ventral)[-c(1,2,3,4)])
plotData_ventral$Significant <- ifelse(plotData_ventral$FDR < 0.05, "Yes", "No")
plotData_ventral$Chr14 <- ifelse(plotData_ventral$genomeID == "chr14", "Chr14", "Other")
plotData_ventral <- arrange(plotData_ventral, desc(PValue))

p_ventral <- ggparcoord(plotData_ventral, groupColumn = 34,
           columns = 1:4,
           alphaLines = 1,
           showPoints = TRUE,
           scale = "globalminmax") +
  xlab("") +
  ylab("logFC") +
  ggtitle("Ventral: RTT vs IC") +
  scale_color_manual(values = c("#69b3a2", "#E8E8E8")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))

plotData_dorsal <- inner_join(topTable_dorsal, featureInfo, by = c("ID" = "miRNA_Index"))
colnames(plotData_dorsal) <- c("D0", "D13", "D40", "D75", colnames(plotData_dorsal)[-c(1,2,3,4)])
plotData_dorsal$Significant <- ifelse(plotData_dorsal$FDR < 0.05, "Yes", "No")
plotData_dorsal$Chr14 <- ifelse(plotData_dorsal$genomeID == "chr14", "Chr14", "Other")
plotData_dorsal <- arrange(plotData_dorsal, desc(PValue))

p_dorsal <- ggparcoord(plotData_dorsal, groupColumn = 34,
                       columns = 1:4,
                       alphaLines = 1,
                       showPoints = TRUE,
                       scale = "globalminmax") +
  xlab("") +
  ylab("logFC") +
  ggtitle("Dorsal: RTT vs IC") +
  scale_color_manual(values = c("#69b3a2", "#E8E8E8")) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))



p_all <- ggarrange(p_ventral, p_dorsal, nrow = 2, ncol = 1, align = "v", common.legend = TRUE, legend = "bottom")
ggsave(p_all, file = "parallelPlot.png", width = 8, height = 6)





plotData_ventral <- inner_join(topTable_ventral, featureInfo, by = c("ID" = "miRNA_Index"))
colnames(plotData_ventral) <- c("D0", "D13", "D40", "D75", colnames(plotData_ventral)[-c(1,2,3,4)])
plotData_ventral$Significant <- ifelse(plotData_ventral$FDR < 0.05, "Yes", "No")
plotData_ventral$Chr14 <- ifelse(plotData_ventral$genomeID == "chr14", "Chr14", "Other")
plotData_ventral <- arrange(plotData_ventral, desc(PValue))

p_ventral <- ggparcoord(plotData_ventral, groupColumn = 33,
                        columns = 1:4,
                        alphaLines = 1,
                        showPoints = TRUE,
                        scale = "globalminmax") +
  xlab("") +
  ylab("logFC") +
  ggtitle("Ventral: RTT vs IC") +
  scale_color_manual(values = rev(c("#69b3a2", "#E8E8E8"))) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))

plotData_dorsal <- inner_join(topTable_dorsal, featureInfo, by = c("ID" = "miRNA_Index"))
colnames(plotData_dorsal) <- c("D0", "D13", "D40", "D75", colnames(plotData_dorsal)[-c(1,2,3,4)])
plotData_dorsal$Significant <- ifelse(plotData_dorsal$FDR < 0.05, "Yes", "No")
plotData_dorsal$Chr14 <- ifelse(plotData_dorsal$genomeID == "chr14", "Chr14", "Other")
plotData_dorsal <- arrange(plotData_dorsal, desc(PValue))

p_dorsal <- ggparcoord(plotData_dorsal, groupColumn = 33,
                       columns = 1:4,
                       alphaLines = 1,
                       showPoints = TRUE,
                       scale = "globalminmax") +
  xlab("") +
  ylab("logFC") +
  ggtitle("Dorsal: RTT vs IC") +
  scale_color_manual(values = rev(c("#69b3a2", "#E8E8E8"))) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))



p_all <- ggarrange(p_ventral, p_dorsal, nrow = 2, ncol = 1, align = "v", common.legend = TRUE, legend = "bottom")
ggsave(p_all, file = "parallelPlot2.png", width = 8, height = 6)


