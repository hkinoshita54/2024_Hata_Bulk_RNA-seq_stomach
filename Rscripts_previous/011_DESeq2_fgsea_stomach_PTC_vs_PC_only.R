# Load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(DESeq2)
library(fgsea)
library(GSVA)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)


# prepare gene sets for gsea ----
collections <- list()
collections$H <- msigdbr(species = "Mus musculus", category = "H")
collections$BIOCARTA <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "BIOCARTA")
collections$BP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
collections$C6 <- msigdbr(species = "Mus musculus", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})


# Load data ----
cts <- read.delim(file = "data/Bulk_stomach_T1PTC_project.txt", header = TRUE, sep = "\t")
rownames(cts) <- cts$Name
cts <- cts[,-1]
cts <- cts[,c(5,6,9,10)]

# DESeq2 ----

### make coldata
coldata <- data.frame(condition = c("PC", "PC", "PTC", "PTC"))
rownames(coldata) <- colnames(cts)
coldata$condition <- factor(coldata$condition, levels = c("PC", "PTC"))

### make DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~condition)

### pre-fileter rows with very few counts
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

### DE analysis
dds <- DESeq(dds)

### PCA plot
# vst = varianceStabilizingTransformation(dds)
vst = vst(dds)
pcaData <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

# PTC vs PC ----

### define contrast
expgroup <- "PTC"
contgroup <- "PC"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
write.table(res,paste0("results/DEG/stomach/PTC_vs_PC_only/", description, ".txt"), sep ="\t", col.names = T, row.names = F)

### volcano plot
# keyvals.colour <- ifelse(res$padj > 0.01 | abs(res$log2FoldChange) < 2, 'gray', 'orange')
# names(keyvals.colour)[keyvals.colour == 'gray'] <- 'other'
# names(keyvals.colour)[keyvals.colour == 'orange'] <- '|logFC|>2 & padj<0.01'
# EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue',
#                 # xlim = c(-10, 10), ylim = c(0, -log10(10e-100)), legendPosition = "right",
#                 # selectLab = c('Lrg1','Cd38'), drawConnectors = TRUE, widthConnectors = 0.5, arrowheads = FALSE,
#                 # pCutoff = 10e-6, FCcutoff = 2, col = c('gray', 'gray', 'gray', 'orange'),
#                 title = paste0(expgroup, " vs ", contgroup), subtitle = NULL, caption = NULL,
#                 colCustom = keyvals.colour,
#                 pointSize = 1.0, labSize = 6.0, cutoffLineWidth = 0.0, colAlpha = 1)

### fgsea

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

#### run fgsea HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/stomach/PTC_vs_PC_only/GSEA_H_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

fgseaRes$pathway <- sub("HALLMARK_", "", fgseaRes$pathway)
fgseaResUp <- fgseaRes[fgseaRes$padj<0.25 & NES>0,]
fgseaResDn <- fgseaRes[fgseaRes$padj<0.25 & NES<0,]
cols <- brewer.pal(3, "Set1")

ggplot(fgseaResUp, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = cols[1]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = "plots/stomach/PTC_vs_PC_only/", width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = "plots/stomach/PTC_vs_PC_only/", width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 

# gene_set = "HALLMARK_HEME_METABOLISM"
# plotEnrichment(collections$H[[gene_set]], rank) + labs(title=gene_set)
