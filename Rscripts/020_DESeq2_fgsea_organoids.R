# Load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(DESeq2)
library(fgsea)
library(GSVA)
# library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(VennDiagram)
library(scales)

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
cts <- read.delim(file = "data/Bulk_organoids_T1PTC_project.txt", header = TRUE, sep = "\t")
rownames(cts) <- cts$Name
cts <- cts[,-1]


# DESeq2 ----

### make coldata
coldata <- data.frame(condition = c("WT", "WT", "PT", "PT", "PC", "PC", "TC", "TC", "PTC", "PTC"))
rownames(coldata) <- colnames(cts)
coldata$condition <- factor(coldata$condition, levels = c("WT", "PT", "PC", "TC", "PTC"))

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
write.table(res,paste0("results/DEG/organoid/", description, ".txt"), sep ="\t", col.names = T, row.names = F)

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

### get DE up
# logFC <- 2
# FDR <- 0.01
# DEG <- res[res$log2FoldChange > logFC & res$padj < FDR,]
# PTCvsPCup <- rownames(DEG)

### get DE up 500
DEG <- res[order(res$stat, decreasing = TRUE),][1:500,]
# DEG <- res[res$log2FoldChange > 1.5,]
# DEG <- DEG[order(DEG$padj),][1:500,]
PTCvsPCup <- rownames(DEG)

### fgsea

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

#### run fgsea HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_H_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

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
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 

# gene_set = "HALLMARK_HEME_METABOLISM"
# plotEnrichment(collections$H[[gene_set]], rank) + labs(title=gene_set)

#### run fgsea C6
fgseaRes <- fgsea(pathways = collections$C6, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_C6_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

#### run fgsea BIOCARTA
fgseaRes <- fgsea(pathways = collections$BIOCARTA, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_BIOCARTA_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

#### run fgsea GOBP
fgseaRes <- fgsea(pathways = collections$BP, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_GOBP_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


# PTC vs WT ----

### define contrast
expgroup <- "PTC"
contgroup <- "WT"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
write.table(res,paste0("results/DEG/organoid/", description, ".txt"), sep ="\t", col.names = T, row.names = F)

### get DE up
# logFC <- 2
# FDR <- 0.01
# DEG <- res[res$log2FoldChange > logFC & res$padj < FDR,]
# PTCvsWTup <- rownames(DEG)

### get DE up 500
DEG <- res[order(res$stat, decreasing = TRUE),][1:500,]
# DEG <- res[res$log2FoldChange > 1.5,]
# DEG <- DEG[order(DEG$padj),][1:500,]
PTCvsWTup <- rownames(DEG)

### fgsea

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

#### run fgsea HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_H_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

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
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 

#### run fgsea C6
fgseaRes <- fgsea(pathways = collections$C6, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_C6_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

#### run fgsea BIOCARTA
fgseaRes <- fgsea(pathways = collections$BIOCARTA, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_BIOCARTA_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

#### run fgsea GOBP
fgseaRes <- fgsea(pathways = collections$BP, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_GOBP_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


# TC vs WT ----

### define contrast
expgroup <- "TC"
contgroup <- "WT"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
write.table(res,paste0("results/DEG/organoid/", description, ".txt"), sep ="\t", col.names = T, row.names = F)

### get DE up
# logFC <- 2
# FDR <- 0.01
# DEG <- res[res$log2FoldChange > logFC & res$padj < FDR,]
# TCvsWTup <- rownames(DEG)

### get DE up 500
DEG <- res[order(res$stat, decreasing = TRUE),][1:500,]
# DEG <- res[res$log2FoldChange > 1.5,]
# DEG <- DEG[order(DEG$padj),][1:500,]
TCvsWTup <- rownames(DEG)

### fgsea

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

#### run fgsea HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_H_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

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
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 

#### run fgsea C6
fgseaRes <- fgsea(pathways = collections$C6, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_C6_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

#### run fgsea BIOCARTA
fgseaRes <- fgsea(pathways = collections$BIOCARTA, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_BIOCARTA_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

#### run fgsea GOBP
fgseaRes <- fgsea(pathways = collections$BP, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_GOBP_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


# PTC vs PT ----

### define contrast
expgroup <- "PTC"
contgroup <- "PT"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
write.table(res,paste0("results/DEG/organoid/", description, ".txt"), sep ="\t", col.names = T, row.names = F)

### fgsea

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

#### run fgsea HALLMARK
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, paste0("results/DEG/organoid/GSEA_H_", description, ".txt"), sep ="\t", col.names = TRUE, row.names = FALSE)

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
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(expgroup, " vs ", contgroup, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = "plots/organoid", width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300)

# GSVA ----

### column annotation
genotype <- c("WT", "PT", "PC", "TC", "PTC")
anno_col <- data.frame(row.names = colnames(cts), organoid = factor(rep(genotype, each = 2), levels = genotype))

### gsva and heatmap from DESeq2 (normalized by vst)
gsvaPar <- gsvaParam(assay(vst), collections$H, kcdf = "Gaussian")
gsva.es <- gsva(gsvaPar, verbose=FALSE)
rownames(gsva.es) <- sub("HALLMARK_", "", rownames(gsva.es))
st_colors <- rep(hue_pal()(5))
names(st_colors) <- genotype
ann_colors = list(organoid = st_colors)
hm <- pheatmap(gsva.es, annotation_col = anno_col, show_colnames = FALSE, fontsize_row = 7.5, cluster_cols = FALSE,
               color=colorRampPalette(c("blue", "white", "red"))(100),
               annotation_colors = ann_colors)
save_pheatmap_pdf <- function(x, filename, width = 6, height = 6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(hm, "plots/organoid/GSVA_DESeq2-vst.pdf")

### gsva and heatmap from TPM (from the CLC output excel file)
X <- read.delim(file = "data/Bulk_organoids_T1PTC_project_CPM.txt", header = TRUE, sep = "\t") %>% column_to_rownames(var = "Name")
colnames(X) <- c("WT1", "WT2", "PT1", "PT2", "PC1", "PC2", "TC1", "TC2", "PTC1", "PTC2")
keep <- rowSums(X > 0) >= 2
X <- X[keep,]
X <- log2(X+1)
gsvaPar <- gsvaParam(as.matrix(X), collections$H, kcdf = "Gaussian")
gsva.es <- gsva(gsvaPar, verbose=FALSE)
rownames(gsva.es) <- sub("HALLMARK_", "", rownames(gsva.es))
hm <- pheatmap(gsva.es, annotation_col = anno_col, show_colnames = FALSE, fontsize_row = 7.5, cluster_cols = FALSE,
               color=colorRampPalette(c("blue", "white", "red"))(100))
save_pheatmap_pdf(hm, "plots/organoid/GSVA_TPM.pdf")


# common DEs ----
commonDE_org <- Reduce(intersect, list(PTCvsPCup,PTCvsWTup,TCvsWTup))
commonDE_org_woPC <- Reduce(intersect, list(PTCvsWTup,TCvsWTup))
# myCol <- brewer.pal(3, "Pastel2")
# venn.diagram(x = list(PTCvsPCup, PTCvsWTup, TCvsWTup), category.names = c("PTC vs PC" , "PTC vs WT" , "TC vs WT"),
#              filename = 'plots/organoid/venn_diagramm.png', output = TRUE, 
#              imagetype="png", height = 480, width = 480, resolution = 300, compression = "lzw",    # output
#              lwd = 2, lty = 'blank', fill = myCol,    # circles
#              cex = .6, fontface = "bold", fontfamily = "sans",    # numbers
#              cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(-27, 27, 135), 
#              cat.dist = c(0.055, 0.055, 0.085), cat.fontfamily = "sans", rotation = 1    # names
#              )

# common DEs among sto and org ----
commonDE_all <- intersect(commonDE_sto, commonDE_org)
commonDE_all_woPC <- intersect(commonDE_sto_woPC, commonDE_org_woPC)
genes <- c("Lrg1", "Cd38")
genes %in% commonDE_all
genes %in% commonDE_all_woPC

myFill <- c(alpha("#1B9E77",0.5), alpha('#D95F02',0.5), alpha('#7570B3',0.5), alpha('#E7298A',0.5))
myCol <- brewer.pal(4, "Dark2")
venn.diagram(x = list(PTCvsWTup, TCvsWTup, PTCvsWTup_sto, TCvsWTup_sto), 
             category.names = c("PTC_o" , "TC_o", "PTC_s", "TC_s"),
             filename = 'plots/organoid/venn_diagramm_2.png', output = TRUE, disable.logging = TRUE,
             imagetype="png", height = 480, width = 480, resolution = 300, compression = "lzw",    # output
             lwd = 2, fill = myFill, col = myCol,   # circles
             cex = 0.5, fontface = "bold", fontfamily = "sans",    # numbers
             cat.cex = 0.5, cat.fontface = "bold", cat.default.pos = "outer", cat.fontfamily = "sans",   # names
             margin = 0.05
)

DE_list <- qpcR:::cbind.na(PTCvsPCup_sto, PTCvsWTup_sto, TCvsWTup_sto, PTCvsPCup, PTCvsWTup,TCvsWTup,
             commonDE_sto, commonDE_sto_woPC, commonDE_org, commonDE_org_woPC, commonDE_all, commonDE_all_woPC)
DE_list[is.na(DE_list)] <- ""
write.table(DE_list, file = "results/DEG/DEG_list.txt", sep = "\t", row.names = FALSE)



# PC vs WT ----

### define contrast
expgroup <- "PC"
contgroup <- "WT"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
write.table(res,paste0("results/DEG/organoid/", description, ".txt"), sep ="\t", col.names = T, row.names = F)