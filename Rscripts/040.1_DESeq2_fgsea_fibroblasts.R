analysis_step <- "040.1_fibro"

# Load packages ----
library(tidyverse)
library(openxlsx2)
library(msigdbr)
library(DESeq2)
library(fgsea)
library(GSVA)
# library(EnhancedVolcano)
library(scico)
library(RColorBrewer)
library(pheatmap)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
deg_path <- file.path(res_path, "DEG")
fs::dir_create(c(plot_path, res_path, deg_path))

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
cts <- read.delim(file = "data/Bulk_fibro_T1PTC_project.txt", header = TRUE, sep = "\t")
rownames(cts) <- cts$Name
cts <- cts[,-1]


# DESeq2 ----

### make coldata
coldata <- data.frame(condition = c("PC", "PC", "PTC", "PTC", "PTC", "aCD38", "aCD38", "aCD38"))
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

condition_col <- c(
  "PC"  = "#009E73",
  "PTC" = "#CC79A7",
  "aCD38" = "#377EB8"
)


ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% var.")) +
  ylab(paste0("PC2: ",percentVar[2],"% var.")) + 
  coord_fixed() +
  theme_bw(base_size = 14) +
  scale_color_manual("conditions", values = condition_col)
ggsave(filename = "pca.png", path = plot_path, width = 3.5, height = 2, dpi = 300) 

# GSVA ----

### column annotation
anno_col <- data.frame(
  row.names = colnames(cts), 
  conditions = factor(
    c("PC", "PC", "PTC", "PTC", "PTC", "aCD38", "aCD38", "aCD38"), 
    levels = c("PC", "PTC", "aCD38")
  )
)

### gsva and heatmap from DESeq2 (normalized by vst)
gsvaPar <- gsvaParam(assay(vst), collections$H, kcdf = "Gaussian")
gsva.es <- gsva(gsvaPar, verbose=FALSE)
rownames(gsva.es) <- sub("HALLMARK_", "", rownames(gsva.es))

names(condition_col) <- c("PC", "PTC", "aCD38")
ann_col = list(conditioh = condition_col)

hm <- pheatmap(gsva.es, annotation_col = anno_col, show_colnames = FALSE, fontsize_row = 7.5, cluster_cols = FALSE,
               color=colorRampPalette(scico(9, palette = "vik"))(100),
               annotation_colors = ann_col)
save_pheatmap_pdf <- function(x, filename, width = 6, height = 6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(hm, file.path(plot_path, "GSVA_DESeq2-vst.pdf"))

# pairwise comoparison ----
resultsNames(dds)

## aCD38 vs PTC ----
test <- "aCD38"
ref <- "PTC"
desc <- paste0(test, "_vs_", ref)
res <- results(dds, contrast = c("condition", test, ref)) %>% data.frame() %>% rownames_to_column(var = "gene_name")
write_xlsx(res, file.path(deg_path, paste0(desc, ".xlsx")), col.names = T, row.names = F)

aCD38_fib_up <- res %>% 
  arrange(desc(stat)) %>% 
  slice_head(n = 500) %>% 
  pull(gene_name)

aCD38_fib_dn <- res %>% 
  arrange((stat)) %>% 
  slice_head(n = 500) %>% 
  pull(gene_name)

save(aCD38_fib_up, aCD38_fib_dn, file = "aCD38_fib_up_dn.RData")

### fgsea
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write_xlsx(fgseaRes, file.path(res_path, paste0("fgsea_H_", desc, ".xlsx")), col_names = TRUE, row_names = FALSE)

fgseaRes$pathway <- sub("HALLMARK_", "", fgseaRes$pathway)
fgseaResUp <- fgseaRes[fgseaRes$padj<0.25 & NES>0,]
fgseaResDn <- fgseaRes[fgseaRes$padj<0.25 & NES<0,]
cols <- brewer.pal(3, "Set1")

ggplot(fgseaResUp, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = cols[1]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(desc, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Up_", desc, ".png"), path = plot_path, width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(desc, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", desc, ".png"), path = plot_path, width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 


## PTC vs PC ----
test <- "PTC"
ref <- "PC"
desc <- paste0(test, "_vs_", ref)
res <- results(dds, contrast = c("condition", test, ref)) %>% data.frame() %>% rownames_to_column(var = "gene_name")
write_xlsx(res, file.path(deg_path, paste0(desc, ".xlsx")), col.names = T, row.names = F)
#
# aCD38_fib_up <- res %>% 
#   arrange(desc(stat)) %>% 
#   slice_head(n = 500) %>% 
#   pull(gene_name)
# aCD38_fib_dn <- res %>% 
#   arrange((stat)) %>% 
#   slice_head(n = 500) %>% 
#   pull(gene_name)
# save(aCD38_fib_up, aCD38_fib_dn, file = "aCD38_fib_up_dn.RData")

### fgsea
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write_xlsx(fgseaRes, file.path(res_path, paste0("fgsea_H_", desc, ".xlsx")), col_names = TRUE, row_names = FALSE)

fgseaRes$pathway <- sub("HALLMARK_", "", fgseaRes$pathway)
fgseaResUp <- fgseaRes[fgseaRes$padj<0.25 & NES>0,]
fgseaResDn <- fgseaRes[fgseaRes$padj<0.25 & NES<0,]
cols <- brewer.pal(3, "Set1")

ggplot(fgseaResUp, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = cols[1]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(desc, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Up_", desc, ".png"), path = plot_path, width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(desc, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", desc, ".png"), path = plot_path, width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300) 
