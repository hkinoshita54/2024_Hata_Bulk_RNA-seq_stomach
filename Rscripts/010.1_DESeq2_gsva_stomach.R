analysis_step <- "010_stomach"

# Load packages ----
library(tidyverse)
library(ggplot2)
library(msigdbr)
library(DESeq2)
# library(fgsea)
library(GSVA)
# library(EnhancedVolcano)
library(scico)
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

yap_list <- readRDS("RDSfiles/yap_list.RDS")
yap_list <- yap_list[4:5]

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
# normalized_cts <- counts(dds, normalized=TRUE)


### PCA plot
# vst = varianceStabilizingTransformation(dds)
vst = vst(dds)
pcaData <- plotPCA(vst, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

genotype_col <- c(
  "WT"  = "#999999",
  "PT"  = "#56B4E9",
  "PC"  = "#009E73",
  "TC"  = "#E69F00",
  "PTC" = "#CC79A7"
)

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% var.")) +
  ylab(paste0("PC2: ",percentVar[2],"% var.")) + 
  # coord_fixed() +
  theme_bw(base_size = 14) +
  scale_color_manual("Genotype", values = genotype_col)
ggsave(filename = "pca.png", path = plot_path, width = 3.5, height = 2, dpi = 300) 

# GSVA ----

### column annotation
genotype <- c("WT", "PT", "PC", "TC", "PTC")
anno_col <- data.frame(row.names = colnames(cts), stomach = factor(rep(genotype, each = 2), levels = genotype))

### gsva and heatmap from DESeq2 (normalized by vst)
gsvaPar <- gsvaParam(assay(vst), collections$H, kcdf = "Gaussian")
gsva.es <- gsva(gsvaPar, verbose=FALSE)
rownames(gsva.es) <- sub("HALLMARK_", "", rownames(gsva.es))

st_colors <- genotype_col
names(st_colors) <- genotype
ann_colors = list(stomach = st_colors)

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
save_pheatmap_pdf(hm, file.path(plot_path, "GSVA_DESeq2-vst.pdf"))

### heatmap of selected gene sets
rownames(gsva.es)
gene_sets_sel <- c(
  "EPITHELIAL_MESENCHYMAL_TRANSITION",
  "ANGIOGENESIS",
  "INFLAMMATORY_RESPONSE",
  "TNFA_SIGNALING_VIA_NFKB",
  "TGF_BETA_SIGNALING",
  "HYPOXIA",
  "INTERFERON_ALPHA_RESPONSE",
  "DNA_REPAIR",
  "FATTY_ACID_METABOLISM",
  "OXIDATIVE_PHOSPHORYLATION"
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
  width = 6, 
  height = 2.5
)

### YAP signatures
gsvaPar_yap <- gsvaParam(assay(vst), yap_list, kcdf = "Gaussian")
gsva.es_yap <- gsva(gsvaPar_yap, verbose=FALSE)

hm <- pheatmap(
  gsva.es_yap, 
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
  file.path(plot_path, "GSVA_yap_DESeq2-vst.pdf"),
  width = 4, 
  height = 1
)


# export DEG ----
get_res <- function(
    test = test,
    ref = ref,
    out_path = out_path
){
  # define contrast
  test <- test
  ref <- ref
  desc <- paste0(test, "_vs_", ref)
  
  # write result to excel
  res <- results(dds, contrast = c("condition", test, ref)) %>% data.frame()
  res <- data.frame(gene_name = rownames(res), res)
  write_xlsx(res, file.path(out_path, paste0(desc, ".xlsx")), col.names = T, row.names = F)
}

get_res(test = "PTC", ref = "PC", out_path = deg_path)
get_res(test = "PTC", ref = "WT", out_path = deg_path)
get_res(test = "PC", ref = "WT", out_path = deg_path)
