# Load packages
library(tidyverse)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rWikiPathways)

# Set working directory
setwd("E:/RTTproject/ExosomeAnalysis/PathwayAnalysis")

# Load table (miRtarbase v9.0)
MTI <- read_excel("hsa_MTI.xlsx")

# Load data (feature info)
homeDir <- "E:/RTTproject/ExosomeAnalysis/"
load(paste0(homeDir,"featureInfo.RData"))
featureInfo <- featureInfo[!duplicated(featureInfo$miR_name),]
featureInfo$CompleteName <- str_remove_all(featureInfo$miR_name, "_.*")

# Filter MTI
MTI_filtered <- MTI[MTI$miRNA %in% featureInfo$CompleteName,]
save(MTI_filtered, file = "MTI_filtered.RData")

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


###############################################################################

# Get Data

###############################################################################

#==============================================================================#
# Early expression
#==============================================================================#
# Get early expressed miRNAs
earlyExpr <- read.table("E:/RTTproject/ExosomeAnalysis/Time Analysis/DownExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% earlyExpr$miRNA,]

#==============================================================================#
# Late expression
#==============================================================================#
# Get late expressed miRNAs
earlyExpr <- read.table("E:/RTTproject/ExosomeAnalysis/Time Analysis/UpExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% earlyExpr$miRNA,]

#==============================================================================#
# Constant expression
#==============================================================================#

# Get constant expressed miRNAs
earlyExpr <- read.table("E:/RTTproject/ExosomeAnalysis/Time Analysis/ConstantExpr.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% earlyExpr$miRNA,]

#==============================================================================#
# RTT vs IC
#==============================================================================#

# Get RTTvsIC miRNAs
earlyExpr <- read.table("E:/RTTproject/ExosomeAnalysis/RTTvsIC/RTTvsIC.txt", col.names = "miRNA")

# Get the associated genes
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% earlyExpr$miRNA,]


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

result_WP_early <- WP@result
result_WP_early_copy <- WP@result[,c(1,2,5)]

# Number of permutations
nPerm <- 10000

# All miRNAs
all_mirnas <- unique(MTI_filtered$miRNA)

# Number of miRNAs to select during each permutation
n_miRNA <- length(earlyExpr$miRNA)

# i = 9431
#set.seed(59807)
set.seed(337872)
#set.seed(222222)
for (i in 4628:nPerm){
  
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
  result_WP_early_copy <- left_join(result_WP_early_copy, results, by = c("ID" = "ID"))
}
save(result_WP_early_copy, file = "results_WP_constant_copy.RData")

# Get permutation results
perm_results <- -log10(result_WP_early_copy[,4:(nPerm+3)])
perm_results[is.na(perm_results)] <- 0

# Count number of times with more enrichment
test <- rowSums(perm_results > -log10(result_WP_early_copy[,3]))/nPerm

# Combine into data.frame
Output <- cbind.data.frame(result_WP_early_copy[,1:3], test)

# Calculate FDR
Output$FDR <- p.adjust(test, method = "fdr")

# Change column names
colnames(Output) <- c("ID", "Description", "pvalue (no perm)", "pvalue (perm)", "FDR")

# Order the rows by pvalue
Output <- arrange(Output, by = `pvalue (perm)`)

# Write output
save(perm_results, file = "perm_results_constant.RData")
write.csv(Output, file = "WP_Constant.csv")




load("perm_results_RTTvsIC.RData")




# Filter GMT
included_pathways <- NULL
uniqueWP <- unique(gmt$wpid)
for (i in 1:length(uniqueWP)) {
  OntTerms <- rWikiPathways::getOntologyTermNames(uniqueWP[i])
  if ((sum(str_detect(OntTerms, "cancer")) == 0) &
      (sum(str_detect(OntTerms, "disease")) == 0)){
    included_pathways <- c(included_pathways, uniqueWP[i])
  }
}

# also include RTT pathways
included_pathways <- c(included_pathways, "WP3584", "WP4312")
gmt <- gmt[gmt$wpid %in% included_pathways,]
#save(gmt, file = "gmt_WP.RData")
load("gmt_WP.RData")








#*****************************************************************************#
# Early expression
#*****************************************************************************#
earlyExpr <- read.table("E:/RTTproject/ExosomeAnalysis/Time Analysis/EarlyExpr.txt", col.names = "miRNA")
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% earlyExpr$miRNA,]

#===============#
# Gene Ontology
#===============#

# Common genes
GOterms_BP <- enrichGO(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
                       universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                       keyType = "ENTREZID",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = Inf,
                       qvalueCutoff = Inf,
                       minGSSize = 10,
                       maxGSSize = 100)

results_BP_early <- GOterms_BP@result
#test <- as.numeric(lapply(str_split(results_BP_early$BgRatio, "/"), `[[`, 1))
#results_BP_early <- results_BP_early[test < 100,]

write.csv(results_BP_early, file = "GO_early.csv")

# Common terms
total <- NULL
for (i in 1:length(earlyExpr$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% earlyExpr$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    GO <- enrichGO(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   keyType = "ENTREZID",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 100)
    
    sig_GO <- GO@result[GO@result$p.adjust < 0.05,]
    #test <- as.numeric(lapply(str_split(sig_GO$BgRatio, "/"), `[[`, 1))
    #sig_GO <- sig_GO[test < 100,]
    
    sig_GO$miRNA <- rep(earlyExpr$miRNA[i],nrow(sig_GO))
    
    total <- rbind.data.frame(total, sig_GO)
  }
  
}

Terms <- names(table(total$ID))[which(table(total$ID) > 1)]

GOHeatmap_early <- ggplot() +
  geom_tile(data = total[total$ID %in% Terms,], 
            aes(y = Description, x = factor(miRNA, levels = earlyExpr$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(earlyExpr$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")
  

ggsave(GOHeatmap_early, file = "GOHeatmap_early.png", width = 8, height = 6)


#===============#
# WikiPathways
#===============#

# common genes
WP <- enricher(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
               universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               TERM2GENE = path2gene,
               TERM2NAME = path2name)

result_WP_early <- WP@result
write.csv(result_WP_early, file = "WP_early.csv")


# common pathways
total <- NULL
for (i in 1:length(earlyExpr$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% earlyExpr$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    WP <- enricher(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   TERM2GENE = path2gene,
                   TERM2NAME = path2name)
    
    sig_WP <- WP@result[WP@result$p.adjust < 0.05,]
    sig_WP$miRNA <- rep(earlyExpr$miRNA[i],nrow(sig_WP))
    
    total <- rbind.data.frame(total, sig_WP)
  }
 
}

Pathways <- names(table(total$ID))[which(table(total$ID) > 2)]

pathHeatmap_early <- ggplot() +
  geom_tile(data = total[total$ID %in% Pathways,], 
            aes(y = Description, x = factor(miRNA, levels = earlyExpr$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(earlyExpr$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

ggsave(pathHeatmap_early, file = "WPHeatmap_early.png", width = 8, height = 6)

#*****************************************************************************#
# Late expression
#*****************************************************************************#

# Get gene set
lateExpr <- read.table("E:/RTTproject/ExosomeAnalysis/Time Analysis/lateExpr.txt", col.names = "miRNA")
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% lateExpr$miRNA,]


#===============#
# Gene Ontology
#===============#

# Common genes
GOterms_BP <- enrichGO(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
                       universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                       keyType = "ENTREZID",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = Inf,
                       qvalueCutoff = Inf,
                       minGSSize = 10,
                       maxGSSize = 100)

results_BP_late <- GOterms_BP@result
#test <- as.numeric(lapply(str_split(results_BP_late$BgRatio, "/"), `[[`, 1))
#results_BP_late <- results_BP_late[test < 100,]


write.csv(results_BP_late, file = "GO_late.csv")


# Common terms
total <- NULL
for (i in 1:length(lateExpr$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% lateExpr$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    GO <- enrichGO(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   keyType = "ENTREZID",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 100)
    
    sig_GO <- GO@result[GO@result$p.adjust < 0.05,]
    #est <- as.numeric(lapply(str_split(sig_GO$BgRatio, "/"), `[[`, 1))
    #sig_GO <- sig_GO[test < 100,]
    
    sig_GO$miRNA <- rep(lateExpr$miRNA[i],nrow(sig_GO))
    
    total <- rbind.data.frame(total, sig_GO)
  }
  
}

Terms <- names(table(total$ID))[which(table(total$ID) > 2)]

GOHeatmap_late <- ggplot() +
  geom_tile(data = total[total$ID %in% Terms,], 
            aes(y = Description, x = factor(miRNA, levels = lateExpr$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(lateExpr$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

ggsave(GOHeatmap_late, file = "GOHeatmap_late.png", width = 8, height = 6)


#===============#
# WikiPathways
#===============#

# Common genes
WP <- enricher(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
               universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               TERM2GENE = path2gene,
               TERM2NAME = path2name)

result_WP_late <- WP@result

write.csv(result_WP_late, file = "WP_late.csv")


# Common pathways
total <- NULL
for (i in 1:length(lateExpr$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% lateExpr$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    WP <- enricher(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   TERM2GENE = path2gene,
                   TERM2NAME = path2name)
    
    sig_WP <- WP@result[WP@result$p.adjust < 0.05,]
    sig_WP$miRNA <- rep(lateExpr$miRNA[i],nrow(sig_WP))
    
    total <- rbind.data.frame(total, sig_WP)
  }
  
}

Pathways <- names(table(total$ID))[which(table(total$ID) > 3)]

pathHeatmap_late <- ggplot() +
  geom_tile(data = total[total$ID %in% Pathways,], 
            aes(y = Description, x = factor(miRNA, levels = lateExpr$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(lateExpr$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

ggsave(pathHeatmap_late, file = "WPHeatmap_late.png", width = 8, height = 6)

#*****************************************************************************#
# Constant expression
#*****************************************************************************#

# Get gene set
constantExpr <- read.table("E:/RTTproject/ExosomeAnalysis/Time Analysis/ConstantExpr.txt", col.names = "miRNA")
geneSet <- MTI_filtered[MTI_filtered$miRNA %in% constantExpr$miRNA,]

# Plot number of targets
table(geneSet$miRNA)
hist(table(MTI_filtered$miRNA))

plotDF_all <- data.frame(nTarget = table(MTI_filtered$miRNA))
plotDF_constant <- data.frame(nTarget = table(geneSet$miRNA))
p <- ggplot() +
  geom_histogram(data = plotDF_all, aes(x = nTarget.Freq), bins = 100) +
  geom_vline(data = plotDF_constant, aes(xintercept = nTarget.Freq), linetype = "dashed",
            color = "red") +
  xlab("Number of target genes") +
  ylab("Count") +
  theme_classic()

ggsave(p, file = "nTargets_constant.png", width = 8, height = 6)
#===============#
# Gene Ontology
#===============#

# Common genes
GOterms_BP <- enrichGO(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
                       universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                       keyType = "ENTREZID",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       minGSSize = 10,
                       maxGSSize = 100)

results_BP_constant <- GOterms_BP@result
#test <- as.numeric(lapply(str_split(results_BP_constant$BgRatio, "/"), `[[`, 1))
#results_BP_constant<- results_BP_constant[test < 100,]


write.csv(results_BP_constant, file = "GO_constant.csv")


# Common terms
total <- NULL
for (i in 1:length(constantExpr$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% constantExpr$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    GO <- enrichGO(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   keyType = "ENTREZID",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 100)
    
    sig_GO <- GO@result[GO@result$p.adjust < 0.05,]
    #test <- as.numeric(lapply(str_split(sig_GO$BgRatio, "/"), `[[`, 1))
    #sig_GO <- sig_GO[test < 100,]
    
    sig_GO$miRNA <- rep(constantExpr$miRNA[i],nrow(sig_GO))
    
    total <- rbind.data.frame(total, sig_GO)
  }
  
}

Terms <- names(table(total$ID))[which(table(total$ID) > 1)]

GOHeatmap_constant <- ggplot() +
  geom_tile(data = total[total$ID %in% Terms,], 
            aes(y = Description, x = factor(miRNA, levels = constantExpr$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(constantExpr$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

ggsave(GOHeatmap_constant, file = "GOHeatmap_constant.png", width = 8, height = 6)


#===============#
# WikiPathways
#===============#

# Common genes
WP <- enricher(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
               universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               TERM2GENE = path2gene,
               TERM2NAME = path2name)

result_WP_constant <- WP@result

write.csv(result_WP_constant, file = "WP_constant.csv")


# Common pathways
total <- NULL
for (i in 1:length(constantExpr$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% constantExpr$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    WP <- enricher(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   TERM2GENE = path2gene,
                   TERM2NAME = path2name)
    
    sig_WP <- WP@result[WP@result$p.adjust < 0.05,]
    sig_WP$miRNA <- rep(constantExpr$miRNA[i],nrow(sig_WP))
    
    total <- rbind.data.frame(total, sig_WP)
  }
  
}

Pathways <- names(table(total$ID))[which(table(total$ID) > 1)]

pathHeatmap_constant <- ggplot() +
  geom_tile(data = total[total$ID %in% Pathways,], 
            aes(y = Description, x = factor(miRNA, levels = constantExpr$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(constantExpr$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

ggsave(pathHeatmap_constant, file = "WPHeatmap_constant.png", width = 8, height = 6)



###############################################################################

# 2. RTT vs IC

###############################################################################

rm(list = ls())

# Load mirRNA-target interaction data
load("MTI_filtered.RData")

#read GMT file
load("gmt_WP.RData")

# link pathway to gene
path2gene <- gmt[,c("wpid", "gene")]
path2name <- gmt[,c("wpid", "name")]

# Load data
miRNAs <- read.table("E:/RTTproject/ExosomeAnalysis/RTTvsIC/RTTvsIC.txt", col.names = "miRNA")

geneSet <- MTI_filtered[MTI_filtered$miRNA %in% miRNAs$miRNA,]


#===============#
# Gene Ontology
#===============#

# Common genes
GOterms_BP <- enrichGO(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
                       universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                       keyType = "ENTREZID",
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       minGSSize = 10,
                       maxGSSize = 100)

results_BP <- GOterms_BP@result
#test <- as.numeric(lapply(str_split(results_BP$BgRatio, "/"), `[[`, 1))
#results_BP <- results_BP[test < 100,]

write.csv(results_BP, file = "GO_RTTvsIC.csv")


# Common terms
total <- NULL
for (i in 1:length(lateExpr$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% lateExpr$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    GO <- enrichGO(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   keyType = "ENTREZID",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   minGSSize = 10,
                   maxGSSize = 100)
    
    sig_GO <- GO@result[GO@result$p.adjust < 0.05,]
    sig_GO$miRNA <- rep(lateExpr$miRNA[i],nrow(sig_GO))
    
    total <- rbind.data.frame(total, sig_GO)
  }
  
}

Terms <- names(table(total$ID))[which(table(total$ID) > 3)]

GOHeatmap_late <- ggplot() +
  geom_tile(data = total[total$ID %in% Terms,], 
            aes(y = Description, x = factor(miRNA, levels = lateExpr$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(lateExpr$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

ggsave(GOHeatmap_late, file = "GOHeatmap_late.png", width = 8, height = 6)


#===============#
# WikiPathways
#===============#

# Common genes
WP <- enricher(gene = as.character(unique(geneSet$`Target Gene (Entrez ID)`)),
               universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               TERM2GENE = path2gene,
               TERM2NAME = path2name)

result_WP <- WP@result

write.csv(result_WP, file = "WP_RTTvsIC.csv")


# Common pathways
total <- NULL
for (i in 1:length(miRNAs$miRNA)){
  geneSet_i <- MTI_filtered[MTI_filtered$miRNA %in% miRNAs$miRNA[i],]
  
  if (nrow(geneSet_i) > 0){
    WP <- enricher(gene = as.character(unique(geneSet_i$`Target Gene (Entrez ID)`)),
                   universe = as.character(unique(MTI_filtered$`Target Gene (Entrez ID)`)),
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 1,
                   TERM2GENE = path2gene,
                   TERM2NAME = path2name)
    
    sig_WP <- WP@result[WP@result$p.adjust < 0.05,]
    sig_WP$miRNA <- rep(miRNAs$miRNA[i],nrow(sig_WP))
    
    total <- rbind.data.frame(total, sig_WP)
  }
  
}

Pathways <- names(table(total$ID))[which(table(total$ID) > 5)]

pathHeatmap <- ggplot() +
  geom_tile(data = total[total$ID %in% Pathways,], 
            aes(y = Description, x = factor(miRNA, levels = miRNAs$miRNA), fill = -log(pvalue))) +
  scale_x_discrete("miRNA", breaks = factor(miRNAs$miRNA), drop = FALSE) +
  scale_fill_viridis_c(option = "magma") +
  labs(fill = "-log(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        axis.title = element_blank(),
        legend.position = "right")

ggsave(pathHeatmap, file = "WPHeatmap_RTTvsIC.png", width = 10, height = 6)


