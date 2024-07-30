####
# generate plot for GeneRank by taiji
# combine it with expression levels and DE analysis
# PTC vs PC (organoids)


# load packages ----
library(biomaRt)
library(org.Mm.eg.db)
library(tidyverse)


# map of human gene symbols to mouse gene symbols ----
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=mart) %>% distinct() %>% as_tibble() 
bm[bm == ""] <- NA
bm <- bm %>% na.omit()


# load data and connect with DE analysis ----
GeneRank <- read.delim(file = "data/Score_Rank_PTCvsPC.txt")    # input is prepared by excel
GeneRank <- inner_join(GeneRank, bm, by = c("TFs"="hsapiens_homolog_associated_gene_name"))    # add mouse gene symbols

####
# do part of 025_DESeq2_commonDE_organoid_logFC>0.58_FDR<0.1.R

GeneRank <- inner_join(GeneRank, res, by = c("external_gene_name"="gene_name"))    # connect to DE analysis, res is from DESeq2
# connect to gene expression level
vsd <- vst(dds, blind = FALSE)    # retrieve transformed noramlized count
mtx <- assay(vsd) %>% as.data.frame() %>% rownames_to_column(var = "gene_name")
GeneRank <- inner_join(GeneRank, mtx, by = c("external_gene_name"="gene_name"))

# plot of GeneRank_logFC vs rank
ggplot(GeneRank, aes(x = rank_PTCvsPC, y = log2FC_PTCvsPC)) +
  geom_point(size = 1, alpha = 1, col = ifelse(GeneRank$log2FoldChange > log2(1.5) & GeneRank$log2FC_PTCvsPC > 0, "red","black")) + 
  geom_text_repel(aes(label = ifelse(log2FoldChange > log2(1.5) & log2FC_PTCvsPC > 0, external_gene_name,"")), colour = "red", size = 3, fontface = "italic") +
  geom_hline(yintercept = 0.0, linetype = "dashed") +
  theme_classic()

# plot of GeneRank_logFC vs DE_logFC
ggplot(GeneRank, aes(x = log2FoldChange, y = log2FC_PTCvsPC)) +
  geom_point(size = 1, alpha = 1, col = ifelse(GeneRank$log2FoldChange > log2(1.5) & GeneRank$log2FC_PTCvsPC > 0, "red","black")) + 
  geom_text_repel(aes(label = ifelse(log2FoldChange > log2(1.5) & log2FC_PTCvsPC > 0, external_gene_name,"")), colour = "red", size = 4, fontface = "italic") +
  geom_hline(yintercept = 0.0, linetype = "dashed") +
  theme_classic()

# plot of GeneRank_logFC vs expression (normalized and transformed by vst)
ggplot(GeneRank, aes(x = (PTC1 + PTC2)/2, y = log2FC_PTCvsPC)) +
  geom_point(size = 1, alpha = 1, col = ifelse(GeneRank$log2FoldChange > log2(1.5) & GeneRank$log2FC_PTCvsPC > 0, "red","black")) + 
  geom_text_repel(aes(label = ifelse(log2FoldChange > log2(1.5) & log2FC_PTCvsPC > 0, external_gene_name,"")), colour = "red", size = 4, fontface = "italic") +
  geom_hline(yintercept = 0.0, linetype = "dashed") +
  theme_classic()
