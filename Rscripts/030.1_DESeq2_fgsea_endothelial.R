analysis_step <- "030_endothelial"

# Load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(DESeq2)
library(GSVA)
library(pheatmap)
library(scico)
library(EnhancedVolcano)

# Make directories ----
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

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
cts <- read.delim(file = "data/Bulk_endothelial.txt", header = TRUE, sep = "\t")
rownames(cts) <- cts$Name
cts <- cts[,-1]


# DESeq2 ----

### make coldata
coldata <- data.frame(condition = c("CD105lo", "CD105lo", "CD105hi", "CD105hi"))
rownames(coldata) <- colnames(cts)
coldata$condition <- factor(coldata$condition, levels = c("CD105lo", "CD105hi"))

### make DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~condition)

### pre-fileter rows with very few counts
smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

### DE analysis
dds <- DESeq(dds)

### PCA plot
vst = vst(dds)
pcaData <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ec_cols <- c("#377EB8", "#E41A1C")
names(ec_cols) <- c("CD105lo", "CD105hi")

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% var.")) +
  ylab(paste0("PC2: ",percentVar[2],"% var.")) + 
  # coord_fixed() +
  scale_color_manual("EC", values = ec_cols) +
  theme_bw(base_size = 14)
  
ggsave(filename = "pca.png", path = plot_path, width = 3.5, height = 2, dpi = 300) 

### define contrast
expgroup <- "CD105hi"
contgroup <- "CD105lo"
description <- paste0(expgroup, "_vs_", contgroup)

### get result
res <- results(dds, contrast = c("condition", expgroup, contgroup)) %>% data.frame()
res <- data.frame(gene_name = rownames(res), res)
openxlsx2::write_xlsx(
  res,
  file.path(res_path, paste0(description, ".xlsx")), 
  col_names = TRUE, 
  row.names = FALSE
)

#### make gene rank
rank <- res %>% dplyr::select(gene_name, stat) %>% na.omit() %>% distinct() %>% 
  group_by(gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

# GSVA ----

### column annotation
genotype <- c("CD105lo", "CD105hi")
anno_col <- data.frame(row.names = colnames(cts), endothelial = factor(rep(genotype, each = 2), levels = genotype))

### gsva and heatmap from DESeq2 (normalized by vst)
gsvaPar <- gsvaParam(assay(vst), collections$H, kcdf = "Gaussian")
gsva.es <- gsva(gsvaPar, verbose=FALSE)
openxlsx2::write_xlsx(gsva.es, file.path(res_path, "gsva.xlsx"), col_names = TRUE, row_names = TRUE)
rownames(gsva.es) <- sub("HALLMARK_", "", rownames(gsva.es))
ann_colors = list(endothelial = ec_cols)

hm <- pheatmap(gsva.es, annotation_col = anno_col, show_colnames = FALSE, fontsize_row = 7.5, cluster_cols = FALSE,
               color=colorRampPalette(scico(9, palette = "vik"))(100),
               annotation_colors = ann_colors)
save_pheatmap_pdf <- function(x, filename, width = 6, height = 6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(
  hm, 
  file.path(plot_path, "GSVA_DESeq2-vst.pdf"),
  width = 5.5, 
  height = 6
)

### heatmap of selected gene sets
rownames(gsva.es)
gene_sets_sel <- c(
  "MYC_TARGETS_V1",
  "MYC_TARGETS_V2",
  "E2F_TARGETS",
  "G2M_CHECKPOINT",
  "TGF_BETA_SIGNALING",
  "APICAL_JUNCTION",
  "APICAL_SURFACE",
  "INFLAMMATORY_RESPONSE",
  "ANGIOGENESIS",
  "EPITHELIAL_MESENCHYMAL_TRANSITION"
)
gsva.es_sel <- gsva.es[gene_sets_sel,]

hm <- pheatmap(
  gsva.es_sel, 
  annotation_col = anno_col, 
  show_colnames = FALSE, 
  fontsize_row = 7.5, 
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  color=colorRampPalette(scico(9, palette = "vik"))(100),
  annotation_colors = ann_colors
)

save_pheatmap_pdf(
  hm, 
  file.path(plot_path, "GSVA_sel_DESeq2-vst.pdf"),
  width = 5.5 , 
  height = 2.5
)

## volcano
genes_to_label <- c(
  "Lrg1",
  "Eng",
  "Aplnr",
  "Plvap",
  "Emcn",
  "Csf3",
  "Selp"
)

EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  selectLab = genes_to_label,
  pCutoff = 0.05,
  FCcutoff = 1,
  xlim = c(-10, 10),
  ylim = c(0, 60),
  pointSize = 2.0,
  labSize = 8.0,
  col = c("grey70", "grey70", "grey70", "#D55E00"),
  colAlpha = 0.8,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  boxedLabels = TRUE,
  legendPosition = "none"
) +
  labs(
    title = NULL,
    caption = NULL,
    x = "log2 fold change",
    y = "-log10 adjusted P value"
  ) 
ggsave(filename = "volcano.png", path = plot_path, width = 5, height = 5, dpi = 300) 
