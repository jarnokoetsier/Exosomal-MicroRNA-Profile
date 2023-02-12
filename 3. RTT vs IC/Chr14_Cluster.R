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
setwd("E:/RTTproject/ExosomeAnalysis/RTTvsIC")

# Set homeDir
homeDir <- "E:/RTTproject/ExosomeAnalysis/"

# Load data
load(paste0(homeDir,"NormalizedExpr.RData"))
load(paste0(homeDir,"featureInfo.RData"))
load(paste0(homeDir,"sampleInfo.RData"))

# Make name of samples
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")
featureInfo <- featureInfo[!duplicated(featureInfo$CompleteName),]

# get chr 14 cluster
cluster_chr14 <- featureInfo[featureInfo$genomeID == "chr14",]
cluster_chr14$start <- as.numeric(cluster_chr14$start)
cluster_chr14 <- cluster_chr14[(cluster_chr14$start < 101100000) &
                                 (cluster_chr14$start > 100800000),]
#cluster_chr14 <- cluster_chr14[(cluster_chr14$start < 101100000) &
  #                               (cluster_chr14$start > 101000000),]
cluster_chr14 <- cluster_chr14[cluster_chr14$strand == "+",]

expr_chr14 <- logcpm[rownames(logcpm) %in% cluster_chr14$miRNA_Index,]
df_chr14 <- gather(as.data.frame(expr_chr14), key = "SampleID", value = "Expr")
df_chr14$mir_id <- rep(rownames(expr_chr14), ncol(expr_chr14))
df_chr14 <- inner_join(df_chr14, featureInfo, by = c("mir_id" = "miRNA_Index"))
df_chr14 <- inner_join(df_chr14, sampleInfo, by = c("SampleID" = "SampleID"))
df_chr14$start <- as.numeric(df_chr14$start)
df_chr14$rep_group <- paste(df_chr14$Group, df_chr14$Tissue, df_chr14$Time, df_chr14$miR_name, sep = "_")

# Get average expression
plotData <- df_chr14 %>%
  group_by(rep_group) %>%
  summarize(avgExpr = mean(Expr),
            start = start,
            Group = Group,
            Tissue = Tissue,
            Time = Time,
            SampleID = SampleID,
            miR_name = miR_name)

plotData$TimeGroup <- paste0(plotData$Group, ": ", plotData$Time)

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

p <- ggarrange(p_dorsal, p_ventral, nrow = 2, common.legend = TRUE, legend = "right")

ggsave(p, file = "exprChr14Cluster.png", height = 6, width = 8)

