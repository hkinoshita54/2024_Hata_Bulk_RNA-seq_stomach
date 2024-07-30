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
collections$KEGG <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
collections$REACTOME <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
collections$BP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
collections$C6 <- msigdbr(species = "Mus musculus", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

# Load data ----
cts <- read.delim(file = "data/Bulk_fibro_selected_T1PTC_project.txt", header = TRUE, sep = "\t") # removed PTC-F9 and aCD38-M4
rownames(cts) <- cts$Name
cts <- cts[,-1]


# DESeq2 ----

### make coldata
coldata <- data.frame(condition = c("PC", "PC", "PTC", "PTC", "aCD38", "aCD38"))
rownames(coldata) <- colnames(cts)
coldata$condition <- factor(coldata$condition, levels = c("PC", "PTC", "aCD38"))

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

### heatmap of sample distances
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(cts)
colnames(sampleDistMatrix) <- colnames(cts)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# PTC vs PC ----

### define contrast
expgroup <- "PTC"
contgroup <- "PC"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
write.table(res,paste0("results/DEG/fibroblasts/", description, ".txt"), sep ="\t", col.names = T, row.names = F)

### volcano plot
# keyvals.colour <- ifelse(res$padj > 0.01 | abs(res$log2FoldChange) < 2, 'gray', 'orange')
# names(keyvals.colour)[keyvals.colour == 'gray'] <- 'other'
# names(keyvals.colour)[keyvals.colour == 'orange'] <- '|logFC|>2 & padj<0.01'
# EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue',
#                 # xlim = c(-10, 10), ylim = c(0, -log10(10e-100)), legendPosition = "right",
#                 selectLab = c('Eng','Angpt2'), drawConnectors = TRUE, widthConnectors = 0.5, arrowheads = FALSE,
#                 # pCutoff = 10e-6, FCcutoff = 2, col = c('gray', 'gray', 'gray', 'orange'),
#                 title = paste0(expgroup, " vs ", contgroup), subtitle = NULL, caption = NULL,
#                 colCustom = keyvals.colour,
#                 pointSize = 1.0, labSize = 6.0, cutoffLineWidth = 0.0, colAlpha = 1)

### fgsea

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

#### run fgsea HALLMARK ----
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_H_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

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
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = "plots/fibroblasts", width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = "plots/fibroblasts", width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 

# gene_set = "HALLMARK_HEME_METABOLISM"
# plotEnrichment(collections$H[[gene_set]], rank) + labs(title=gene_set)

#### run fgsea KEGG ----
fgseaRes <- fgsea(pathways = collections$KEGG, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_KEGG_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$KEGG, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_KEGG_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


#### run fgsea REACTOME ----
fgseaRes <- fgsea(pathways = collections$REACTOME, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_REACTOME_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$REACTOME, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_REACTOME_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


#### run fgsea GO:BP ----
fgseaRes <- fgsea(pathways = collections$BP, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_BP_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$BP, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_BP_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


#### run fgsea C6
fgseaRes <- fgsea(pathways = collections$C6, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_C6_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$C6, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_C6_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)



# aCD38 vs PTC ----

### define contrast
expgroup <- "aCD38"
contgroup <- "PTC"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
write.table(res,paste0("results/DEG/fibroblasts/", description, ".txt"), sep ="\t", col.names = T, row.names = F)


### fgsea

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

#### run fgsea HALLMARK ----
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_H_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

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
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = "plots/fibroblasts", width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = "plots/fibroblasts", width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 

# gene_set = "HALLMARK_HEME_METABOLISM"
# plotEnrichment(collections$H[[gene_set]], rank) + labs(title=gene_set)

#### run fgsea KEGG ----
fgseaRes <- fgsea(pathways = collections$KEGG, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_KEGG_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$KEGG, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_KEGG_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


#### run fgsea REACTOME ----
fgseaRes <- fgsea(pathways = collections$REACTOME, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_REACTOME_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$REACTOME, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_REACTOME_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


#### run fgsea GO:BP ----
fgseaRes <- fgsea(pathways = collections$BP, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_BP_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$BP, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_BP_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


#### run fgsea C6
fgseaRes <- fgsea(pathways = collections$C6, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/fibroblasts/GSEA_C6_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
# collapse pathways
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], collections$C6, rank)
fgseaResCollapse <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
write.table(fgseaResCollapse, paste0("results/DEG/fibroblasts/GSEA_C6_Collapse_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)



# GSVA ----

### column annotation
type <- c("PC", "PTC", "aCD38")
anno_col <- data.frame(row.names = colnames(cts), type = factor(c("PC", "PC", "PTC", "PTC", "aCD38", "aCD38"), levels = type))

### gsva and heatmap from DESeq2 (normalized by vst)
gsvaPar <- gsvaParam(assay(vst), collections$H, kcdf = "Gaussian")
gsva.es <- gsva(gsvaPar, verbose=FALSE)
rownames(gsva.es) <- sub("HALLMARK_", "", rownames(gsva.es))
hm <- pheatmap(gsva.es, annotation_col = anno_col, show_colnames = FALSE, fontsize_row = 7.5, cluster_cols = FALSE,
               color=colorRampPalette(c("blue", "white", "red"))(100))
save_pheatmap_pdf <- function(x, filename, width = 6, height = 6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(hm, "plots/fibroblasts/GSVA_DESeq2-vst_selected.pdf")
