# test GCclassifier for mouse RNA-seq data
# https://github.com/Ronlee12355/GCclassifier
analysis_step <- "206_GCclassifier_TPM_log2"

# Load packages ----
library(tidyverse)
library(GCclassifier)
library(biomaRt)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Load data ----
cts <- read.delim(file = "data/TPM_T1PTC.txt", header = TRUE, sep = "\t")
rownames(cts) <- cts$Name
ncts <- cts[,-1]
ncts <- log2(ncts + 1)

# create a table to map mouse gene names to human symbol ----
mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
bm <- getBM(attributes=c("external_gene_name", "hsapiens_homolog_associated_gene_name"), mart=mart) %>%
  distinct() %>%
  as_tibble() 
bm[bm == ""] <- NA
bm <- bm %>% na.omit()

# convert mouse gene symbol to human ----
ncts <- rownames_to_column(ncts, var = "symbol")
ncts <- inner_join(ncts, bm, by = c("symbol"="external_gene_name"))
ncts <- distinct(ncts, hsapiens_homolog_associated_gene_name, .keep_all = TRUE)
ncts <- column_to_rownames(ncts, var = "hsapiens_homolog_associated_gene_name")
ncts <- ncts[, -1]

# GCclassifier
emp.res <- classifyGC(
  Expr = ncts, 
  method = 'EMP',
  idType = 'SYMBOL', 
  useMinPosterior = F 
)
openxlsx2::write_xlsx(emp.res, file.path(res_path, paste0("emp.res", ".xlsx")))

acrg.res <- classifyGC(
  Expr = ncts, 
  method = 'ACRG', 
  idType = 'SYMBOL', 
)
openxlsx2::write_xlsx(acrg.res, file.path(res_path, paste0("acrg.res", ".xlsx")))

tcga.res <- classifyGC(
  Expr = ncts, 
  method = 'TCGA', 
  idType = 'SYMBOL', 
)
openxlsx2::write_xlsx(tcga.res, file.path(res_path, paste0("tcga.res", ".xlsx")))
