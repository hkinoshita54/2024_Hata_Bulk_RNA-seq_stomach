# Load packages ----
library(tidyverse)
library(readxl)
library(tximport)
library(DESeq2)
library(ggplot2)
library(fgsea)
library(RColorBrewer)

# Make directories ----
plot_path <- file.path("plot")
res_path <- file.path("result")
fs::dir_create(plot_path)
fs::dir_create(res_path)

# # Prepare gene sets for gsea ----
# 
# library(biomaRt)
# library(org.Mm.eg.db)
# library(msigdbr)
# 
# collections <- list()
# collections$H <- msigdbr(species = "Mus musculus", category = "H")
# collections$KEGG <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
# collections$REACTOME <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
# collections$BP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
# collections$C6 <- msigdbr(species = "Mus musculus", category = "C6")
# collections <- lapply(collections, function(x) {
#   out <- split(x = x$gene_symbol, f = x$gs_name)
# })
# 
# ensembl <- useEnsembl(biomart = "ensembl", version = 100)
# ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
# anno_mm <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
#                  uniqueRows = TRUE,
#                  mart = ensembl)
# 
# save(anno_mm, collections, file = file.path("RDSfiles", "anno_geneset.RData"))
# ----

# Read data ----
sample_table <- read_excel(file.path("data", "sample_table.xlsx")) %>% as.data.frame
sample_name <- sample_table$sample_name
exp_group <- sample_table$exp_group
rownames(sample_table) <- sample_name

# DESeq2 from tximport
sample_files <- file.path("data", "quants", sample_name, "quant.sf")
names(sample_files) <- sample_name
tx2gene <- read_csv("/Users/gut/ref_genome/mm10/gene_map_mm.csv", col_names = c("ensmustid", "ensmusgid"))
cts = tximport(files = sample_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
dds <- DESeqDataSetFromTximport(txi = cts, colData = sample_table, design = ~exp_group)

# Prefilter genes
smallestGroupSize <- min(table(exp_group))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds$exp_group <- relevel(dds$exp_group, ref = exp_group[1])

# Do DE analysis and save the object
dds <- DESeq(dds)
saveRDS(dds, file = file.path("RDSfiles", "dds.RDS"))

### PCA plot
vsd = vst(dds)
pcaData <- plotPCA(vsd, intgroup = "exp_group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = exp_group, label = sample_name)) +
  geom_text(aes(label = name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()
ggsave(file.path("plot", "PCA.png"), height = 3, width = 5, dpi = 150)

# Pairwise comparison ----
# dds <- readRDS(file.path("RDSfiles", "dds.RDS"))
load(file.path("RDSfiles", "anno_geneset.RData"))
exp_level <- unique(dds$exp_group) %>% as.character()

# Pairwise exp_group 2 vs 1 ----
test_group <- exp_level[2]
ref_group <- exp_level[1]
source(file.path("Rscripts", "050_1_pairwise_comparison.R"))

# Pairwise exp_group 3 vs 1 ----
test_group <- exp_level[3]
ref_group <- exp_level[1]
source(file.path("Rscripts", "050_1_pairwise_comparison.R"))