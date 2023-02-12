
################################################################################

# Read and format data

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(gridExtra)
library(grid)
library(heatmaply)
library(readxl)

# Set home directory
homeDir <- "E:/RTTproject/ExosomeAnalysis/"

# Set working directory
setwd(paste0(homeDir,"ProteinMarkers"))

# Load data
load(paste0(homeDir,"sampleInfo.RData"))

# Read data
pxData_exo <- read.delim("E:/RTTproject/OriginalData/Proteomics/DEP analysis/Result files/20210308_exo_total_scaling_norm.txt")

# Filter data
pxData_exo <- pxData_exo[!pxData_exo$Contaminant,]
pxData_exo <- pxData_exo[pxData_exo$Score.Sequest.HT>10,]
rownames(pxData_exo) <- pxData_exo[,1]

#Remove irrelevant columns
pxData_exo <- pxData_exo[,-c(1:6, 49:59)]

#Remove empty rows
pxData_exo <- pxData_exo[!apply(pxData_exo == "", 1, all),]

#Make dataframe numeric
pxData_exo <- mutate_all(pxData_exo, na_if,"")

for (i in 1:ncol(pxData_exo)){
  pxData_exo[,i] <- str_replace_all(pxData_exo[,i],",", ".")
  pxData_exo[,i] <- as.numeric(pxData_exo[,i])
}

#Make sample IDs correct
samples <- str_remove_all(colnames(pxData_exo), "\\.")
samples <- str_remove(samples, "F..")
samples1 <- rep(0, length(samples))
for (i in 1:length(samples)){
  
  #IC or RTT
  if (str_detect(samples[i],"IC")){
    a <- "IC_"
  } else{
    a <- ""
  }
  
  a <- paste0(a, "MeCP2_R255X_")
  
  #Day
  if (str_detect(samples[i],"D0")){
    a <- paste0(a, "D0")
  } else if (str_detect(samples[i], "D13")){
    a <- paste0(a, "D13")
  } else if (str_detect(samples[i], "D40")){
    a <- paste0(a, "D40")
  } else if (str_detect(samples[i], "D75")){
    a <- paste0(a, "D75")
  }
  
  #Tissue
  if (str_detect(samples[i],"Dorsal")){
    a <- paste0(a, "_Dorsal")
  } else if (str_detect(samples[i], "Ventral")){
    a <- paste0(a, "_Ventral")
  }
  
  #Replicate
  if (str_detect(samples[i],"n1")){
    a <- paste0(a, "_1")
  } else if (str_detect(samples[i], "n2")){
    a <- paste0(a, "_2")
  } else if (str_detect(samples[i], "n3")){
    a <- paste0(a, "_3")
  }
  samples1[i] <- a
}

# Check if all samples are in the data
all(samples1 %in% sampleInfo$SampleID)

# Set column names
colnames(pxData_exo) <- samples1

# log transform data
pxData_exo <- log2(pxData_exo)

# Save data
save(pxData_exo, file = "pxData_exo.RData")


# Read the top 100 exosomal markers from ExoCarta
markers <- read_excel("Exo_markers.xlsx")

# Make annotations: gene symbol to uniprot
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

annotations <- getBM(attributes=c("hgnc_symbol", "protein_id"), 
                     filters = 'hgnc_symbol',
                     values = markers$`Gene Symbol`,
                     mart = ensembl)

write.csv(annotations, file = "annotation.csv")

################################################################################

# Density plot

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)

# Set home directory
homeDir <- "E:/RTTproject/ExosomeAnalysis/"

# Set working directory
setwd(paste0(homeDir,"ProteinMarkers"))

# Load data
load(paste0(homeDir,"sampleInfo.RData"))
load("pxData_exo.RData")

# Gather normalized data
exprPlot <- gather(as.data.frame(pxData_exo), key = "Sample", value = "Expr")

# Combine sampleInfo object with expression data
exprPlot <- inner_join(exprPlot, sampleInfo, by = c("Sample" = "SampleID"))

# Make density plot
exprDensityplot <- ggplot(exprPlot, aes(x = Expr)) +
  geom_density(fill = "#13005A", alpha = 0.8, color = "black", size = 0.8) +
  xlab(expression(log[2]~"intensity")) +
  ylab("Density") +
  theme_minimal() +
  theme(legend.position = "none") 

ggsave(exprDensityplot, file = "exo_proteinExpr_density.png", height = 1, width = 8)



###############################################################################

# Heatmap

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(biomaRt)

# Set homeDir
homeDir <- "E:/RTTproject/ExosomeAnalysis/"

# Set working directory
setwd(paste0(homeDir,"ProteinMarkers"))

# Load sample information data
load(paste0(homeDir,"sampleInfo.RData"))

# Load protein expression data
load("pxData_exo.RData")

# Read the top 100 exosomal markers from ExoCarta
markers <- read_excel("Exo_markers.xlsx")

# Read annotation data (to convert UniProt ID to Protein Name)
annotation <- unique(read.csv("annotation.csv"))

# Remove empty cells
annotation <- annotation[(annotation$UniProtKB.Swiss.Prot.ID != ""),]

# Add missing markers manually in annotation
markers[!(markers$`Gene Symbol` %in% annotation$HGNC.symbol),]
missing <- data.frame(Name = c("HIST1H4A", "HIST2H4A", "HIST1H4B"),
                     ID = c("P62805", "P62805", "P62805"))
colnames(missing) <- colnames(annotation)
annotation <- rbind.data.frame(annotation,missing)

# Get Protein expression
pxData_proteins <- pxData_exo[str_remove(rownames(pxData_exo), "-.") %in% annotation$UniProtKB.Swiss.Prot.ID,]
plotDF <- gather(pxData_proteins)
plotDF$Protein <- rep(rownames(pxData_proteins), ncol(pxData_proteins))
plotDF$ProteinName <- str_remove(plotDF$Protein, "-.")

# Add annotation to protein expression data
plotDF <- inner_join(plotDF, annotation, by = c("ProteinName" = "UniProtKB.Swiss.Prot.ID"))

# Add sample information to protein expression data
plotDF <- inner_join(plotDF, sampleInfo, by = c("key" = "SampleID"))


#*****************************************************************************#
# IC
#*****************************************************************************#

# Filter for IC samples
plotDF_IC <- plotDF[plotDF$Group == "IC",]

# Fix order of gene symbols
plotDF_IC$HGNC.symbol <- factor(plotDF_IC$HGNC.symbol,
                                levels = rev(markers$`Gene Symbol`))

# Correct sample labels
sampleInfo_fil <- sampleInfo[sampleInfo$Group == "IC",]
sampleInfo_fil$Tissue[sampleInfo_fil$Tissue == "Cell"] <- "iPSC"
plotDF_IC$Tissue[plotDF_IC$Tissue == "Cell"] <- "iPSC"
plotDF_IC$key <- factor(plotDF_IC$key,
                        levels = c(rev(sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "Dorsal"]),
                                   sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "iPSC"],
                                   sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "Ventral"]))


# Main plot: heatmap
main <- ggplot(plotDF_IC) +
  geom_tile(aes(x = key, y = HGNC.symbol, fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  labs(fill = expression(log[2] ~ "intensity")) +
  scale_fill_viridis_c(limits = c(15, 30), oob = scales::squish) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 8),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Col side color (Time and tissue)
colSideColor<- unique(data.frame(
  sample = factor(str_remove_all(plotDF_IC$key, "MeCP2_R255X_"),
                  levels = str_remove_all(levels(plotDF_IC$key), "MeCP2_R255X_")),
  time = plotDF_IC$Time,
  tissue = plotDF_IC$Tissue))


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


colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[str_detect(colSideColor$sample, "2") & 
                                  (str_detect(colSideColor$sample, "D40") |
                                     str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  ggtitle("IC-MeCP2:R255X (IC)") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15)) +
  scale_fill_brewer(palette = "Dark2")

# Combine plots
markerPlot <- ggpubr::ggarrange(colSideColorPlot_tissue,
                                main,
                                colSideColorPlot_time,
                                nrow = 3,
                                ncol = 1,
                                heights = c(0.8,9,0.5),
                                align = "v",
                                common.legend = FALSE)

# save figure
ggsave(markerPlot, file = "markerPlot_IC.png", width = 6, height = 9)

#*****************************************************************************#
# RTT
#*****************************************************************************#

# Filter for RTT samples
plotDF_RTT <- plotDF[plotDF$Group == "RTT",]

# Fix order of gene symbols
plotDF_RTT$HGNC.symbol <- factor(plotDF_RTT$HGNC.symbol,
                                levels = rev(markers$`Gene Symbol`))

# Correct sample labels
sampleInfo_fil <- sampleInfo[sampleInfo$Group == "RTT",]
sampleInfo_fil$Tissue[sampleInfo_fil$Tissue == "Cell"] <- "iPSC"
plotDF_RTT$Tissue[plotDF_RTT$Tissue == "Cell"] <- "iPSC"
plotDF_RTT$key <- factor(plotDF_RTT$key,
                        levels = c(rev(sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "Dorsal"]),
                                   sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "iPSC"],
                                   sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "Ventral"]))

# main plot: heatmap
main <- ggplot(plotDF_RTT) +
  geom_tile(aes(x = key, y = HGNC.symbol, fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  labs(fill = expression(log[2] ~ "intensity")) +
  scale_fill_viridis_c(limits = c(15, 30), oob = scales::squish) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 8),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Col side color (Time and tissue)
colSideColor<- unique(data.frame(
  sample = factor(str_remove_all(plotDF_RTT$key, "MeCP2_R255X_"),
                  levels = str_remove_all(levels(plotDF_RTT$key), "MeCP2_R255X_")),
  time = plotDF_RTT$Time,
  tissue = plotDF_RTT$Tissue))


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


colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[str_detect(colSideColor$sample, "2") & 
                                  (str_detect(colSideColor$sample, "D40") |
                                     str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  ggtitle("MeCP2:R255X (RTT)") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 15)) +
  scale_fill_brewer(palette = "Dark2")


# Combine plots
markerPlot <- ggpubr::ggarrange(colSideColorPlot_tissue,
                                main,
                                colSideColorPlot_time,
                                nrow = 3,
                                ncol = 1,
                                heights = c(0.8,9,0.5),
                                align = "v",
                                common.legend = FALSE)

# Save figure
ggsave(markerPlot, file = "markerPlot_RTT.png", width = 6, height = 9)



#*****************************************************************************#
# Make legend
#*****************************************************************************#
legendPlot <- ggplot(plotDF_IC) +
  geom_tile(aes(x = key, y = HGNC.symbol, fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  labs(fill = expression(log[2] ~ "intensity")) +
  scale_fill_viridis_c(limits = c(15, 30), oob = scales::squish) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 8),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank())
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)













mm <- markers$`Gene Symbol`[!(markers$`Gene Symbol` %in% plotDF_IC$HGNC.symbol)]
for (i in mm){
  test <- plotDF_IC
  test[,2:5] <- NA
  test <- unique(test)
  test$HGNC.symbol <- rep(i, nrow(test))
  
  plotDF_IC <- rbind.data.frame(plotDF_IC, test)
}



plotDF_IC$HGNC.symbol <- factor(plotDF_IC$HGNC.symbol,
                                levels = rev(markers$`Gene Symbol`))

sampleInfo_fil <- sampleInfo[sampleInfo$Group == "RTT",]
sampleInfo_fil$Tissue[sampleInfo_fil$Tissue == "Cell"] <- "iPSC"
plotDF_IC$Tissue[plotDF_IC$Tissue == "Cell"] <- "iPSC"
plotDF_IC$key <- factor(plotDF_IC$key,
                            levels = c(rev(sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "Dorsal"]),
                                       sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "iPSC"],
                                       sampleInfo_fil$SampleID[sampleInfo_fil$Tissue == "Ventral"]))


main <- ggplot(plotDF_IC) +
  geom_tile(aes(x = key, y = HGNC.symbol, fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  labs(fill = expression(log[2] ~ "intensity")) +
  scale_fill_viridis_c(limits = c(15, 30), oob = scales::squish) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 8),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Col side color
colSideColor<- unique(data.frame(
  sample = factor(str_remove_all(plotDF_IC$key, "MeCP2_R255X_"),
                  levels = str_remove_all(levels(plotDF_IC$key), "MeCP2_R255X_")),
  time = plotDF_IC$Time,
  tissue = plotDF_IC$Tissue))


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


colSideColorPlot_tissue <- ggplot(data = colSideColor) +
  geom_tile(aes(x = sample, y = "label", fill = tissue)) +
  geom_text(data = colSideColor[str_detect(colSideColor$sample, "2") & 
                                  (str_detect(colSideColor$sample, "D40") |
                                     str_detect(colSideColor$sample, "D0")),],
            aes(x = sample, y = "label", label = tissue)) +
  facet_grid(.~tissue, scales = "free", space = "free") +
  ggtitle("Rett Syndrome") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 20)) +
  scale_fill_brewer(palette = "Dark2")


markerPlot <- ggpubr::ggarrange(colSideColorPlot_tissue,
                       main,
                       colSideColorPlot_time,
                       nrow = 3,
                       ncol = 1,
                       heights = c(0.5,9,0.5),
                       align = "v",
                       common.legend = FALSE)

ggsave(markerPlot, file = "markerPlot_RTT.png", width = 6, height = 14)



legendPlot <- ggplot(plotDF_IC) +
  geom_tile(aes(x = key, y = HGNC.symbol, fill = value)) +
  xlab(NULL) +
  ylab(NULL) +
  facet_grid(.~Tissue, scales = "free", space = "free") +
  labs(fill = expression(log[2] ~ "intensity")) +
  scale_fill_viridis_c(limits = c(15, 30), oob = scales::squish) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 8),
        axis.ticks.y = element_line(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank())
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)



