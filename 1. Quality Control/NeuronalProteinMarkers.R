
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(gridExtra)
library(grid)
library(heatmaply)
library(readxl)

# Set working directory
setwd("D:/RTTproject/ExosomeAnalysis")

# Get data directory
dataDir <- "D:/RTTproject/ExosomeAnalysis/"

# Get output directory
outputDir <- "D:/RTTproject/ExosomeAnalysis/Quality Control/"

# Load exosomal protein expression data
load(paste0(dataDir,"pxData_exo.RData"))

# Load sample information
load(paste0(dataDir,"sampleInfo.RData"))


# Neural markers
markers <- c("L1CAM", "NCAM1", "PAX6", "GSX2",
             "CD63", "GRIA1","GRIA2", "GRIA3", "GRIA4",
             "OLIG2", "EFNA5", "HES1", "HES5",
             "GFAP", "GLAST", "SOX2", "FABP7",
             "TNC", "NES", "VIME", "CADH2",
             "MAP2", "RBFOX3/NeuN", "TBR1", "NEUROD1",
             "NFL", "NFM", "NFH", "TAU",
             "THY1", "UCHL1", "TUBB3", "GAD1",
             "GAD2")
marker_ids <- c("P32004", "P13591", "P26367", "Q9BZM3", 
                "P08962", "P42261", "P42262", "P42263", "P48058",
                "Q13516", "P52803", "Q14469", "Q5TA89",
                "P14136", "P43003", "P48431", "O15540",
                "P24821", "P48681", "P08670", "P19022",
                "P11137", "A6NFN3", "Q16650", "Q13562",
                "P07196", "P07197", "P12036", "P10636",
                "P04216", "P09936", "Q13509", "Q99259",
                "Q05329")
sum(str_remove(rownames(pxData_exo),"-.*") %in% marker_ids)

markers[marker_ids %in% str_remove(rownames(pxData_exo),"-.*")]

#NCAM1, vimentin, N-cadherin, THY1, UCHL1, TUBB3

#https://www.cellsignal.com/pathways/neuronal-and-glial-cell-markers