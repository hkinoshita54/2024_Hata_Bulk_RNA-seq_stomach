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


# PTC vs PC ----

### define contrast
expgroup <- "PTC"
contgroup <- "PC"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)

### get DE up with logFC & FDR
logFC <- log2(1.5)
FDR <- 0.1
DEG <- res[res$log2FoldChange > logFC & res$padj < FDR,]
PTCvsPCup <- rownames(DEG)
# 
# ### get DE top n by stat
# topn <- 300
# DEG <- res[order(res$stat, decreasing = TRUE),][1:topn,]
# 
# ### get DE up with logFC & top n
# logFC <- 2
# topn <- 300
# DEG <- res[res$log2FoldChange > logFC,]
# DEG <- DEG[order(DEG$padj),][1:topn,]
# PTCvsPCup <- rownames(DEG)


# PTC vs WT ----

### define contrast
expgroup <- "PTC"
contgroup <- "WT"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)

### get DE up with logFC & FDR
logFC <- log2(1.5)
FDR <- 0.1
DEG <- res[res$log2FoldChange > logFC & res$padj < FDR,]
PTCvsWTup <- rownames(DEG)


# TC vs WT ----

### define contrast
expgroup <- "TC"
contgroup <- "WT"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)

### get DE up with logFC & FDR
logFC <- log2(1.5)
FDR <- 0.1
DEG <- res[res$log2FoldChange > logFC & res$padj < FDR,]
TCvsWTup <- rownames(DEG)


# common DEs ----
commonDE_org <- Reduce(intersect, list(PTCvsPCup, PTCvsWTup,TCvsWTup))
genes <- c("Lrg1", "Cd38")
genes %in% commonDE_org

# common DEs among sto and org ----
commonDE_all <- intersect(commonDE_sto, commonDE_org)
commonDE_all
genes <- c("Lrg1", "Cd38")
genes %in% commonDE_all

myFill <- c(alpha("#1B9E77",0.5), alpha('#D95F02',0.5), alpha('#7570B3',0.5))
myCol <- brewer.pal(3, "Dark2")

venn.diagram(x = list(PTCvsPCup, PTCvsWTup, TCvsWTup), category.names = c("PTC vs PC" , "PTC vs WT" , "TC vs WT"),
             filename = 'plots/organoid/venn_diagramm_org_FC1.5_FDR0.1.png', output = TRUE, disable.logging = TRUE,
             imagetype="png", height = 480, width = 480, resolution = 300, compression = "lzw", margin = 0.05,    # output
             lwd = 2, col = myCol, fill = myFill,    # circles
             cex = .6, fontface = "bold", fontfamily = "sans",    # numbers
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(-27, 27, 135),
             cat.dist = c(0.08, 0.08, 0.11), cat.fontfamily = "sans", rotation = 1    # names
             )

DE_list <- qpcR:::cbind.na(PTCvsPCup_sto, PTCvsWTup_sto, TCvsWTup_sto, PTCvsPCup, PTCvsWTup,TCvsWTup,
                           commonDE_sto, commonDE_org, commonDE_all)
DE_list[is.na(DE_list)] <- ""
write.table(DE_list, file = "results/DEG/DEG_list_FC1.5_FDR0.1.txt", sep = "\t", row.names = FALSE)
