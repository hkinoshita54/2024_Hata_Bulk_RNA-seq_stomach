# Load packages ----
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(VennDiagram)


# Load data ----
cts <- read.delim(file = "data/Bulk_stomach_T1PTC_project.txt", header = TRUE, sep = "\t")
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
PTCvsPCup_sto <- rownames(DEG)

# ### get DE top n by stat
# topn <- 300
# DEG <- res[order(res$stat, decreasing = TRUE),][1:topn,]
# 
# ### get DE up with logFC & top n
# logFC <- 2
# topn <- 300
# 
# DEG <- res[res$log2FoldChange > 2,]
# DEG <- DEG[order(DEG$padj),][1:300,]
# PTCvsPCup_sto <- rownames(DEG)


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
PTCvsWTup_sto <- rownames(DEG)


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
TCvsWTup_sto <- rownames(DEG)


# common DEs ----
commonDE_sto <- Reduce(intersect, list(PTCvsPCup_sto, PTCvsWTup_sto,TCvsWTup_sto))
genes <- c("Lrg1", "Cd38")
genes %in% commonDE_sto

myFill <- c(alpha("#1B9E77",0.5), alpha('#D95F02',0.5), alpha('#7570B3',0.5))
myCol <- brewer.pal(3, "Dark2")

venn.diagram(x = list(PTCvsWTup_sto, TCvsWTup_sto, PTCvsPCup_sto), category.names = c("PTC vs WT" , "TC vs WT", "PTC vs PC"),
             filename = 'plots/organoid/venn_diagramm_sto_FC1.5_FDR0.1.png', output = TRUE, disable.logging = TRUE,
             imagetype="png", height = 480, width = 480, resolution = 300, compression = "lzw", margin = 0.05,    # output
             lwd = 2, col = myCol, fill = myFill,    # circles
             cex = .6, fontface = "bold", fontfamily = "sans",    # numbers
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer", cat.pos = c(-27, 27, 135),
             cat.dist = c(0.08, 0.08, 0.11), cat.fontfamily = "sans", rotation = 1    # names
)
