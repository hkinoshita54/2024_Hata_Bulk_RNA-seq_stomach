## Load packages ----
library(biomaRt)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(fgsea)

## Load data ----
cts <- read.delim(file = "data/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE, sep = "\t")
rownames(cts) <- cts$gene_name
cts <- cts[,-c(1:2)]

## DESeq2 ----

### make coldata from sra_table
coldata <- data.frame(condition = c("PC", "PC", "PTC", "PTC"))
rownames(coldata) <- colnames(cts)
coldata$condition <- factor(coldata$condition, levels = c("PC", "PTC"))

### make sure cts and coldata are in the same order
# all(rownames(coldata) %in% colnames(cts))
# all(rownames(coldata) == colnames(cts))
# cts <- cts[, rownames(coldata)]
# all(rownames(coldata) == colnames(cts))

### make DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~condition)

#### pre-fileter rows with very few counts
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

#### set factor levels of condition
# dds$condition <- factor(dds$condition, levels = c("WT","PC", "PTC"))

### DE analysis
dds <- DESeq(dds)

### PCA plot
vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

## PTC vs PC ----

### get result, volcano plot
res <- results(dds, contrast = c("condition", "PTC", "PC")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(res,"results/DEG/organoid_PTC_vs_PC_star-salmon.txt", sep ="\t", col.names = T, row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.05, gene_name,"")), colour = "red", size = 3)

### GSEA by fgsea

#### make gene rank
rank <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
rank <- deframe(rank)

#### prepare gene sets
collections <- list()
# collections$GO_BP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
# collections$CGP <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Mus musculus", category = "H")
# collections$C6 <- msigdbr(species = "Mus musculus", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})

#### run fgsea
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes,"results/DEG/organoid_GSEA_H_PTC_vs_PC_salmon.txt", sep ="\t", col.names = TRUE, row.names = FALSE)
