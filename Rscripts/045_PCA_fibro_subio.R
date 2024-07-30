# Load packages ----
library(tidyverse)
library(ggplot2)
library(ggfortify)


cts <- read.delim(file = "data/Bulk_fibro_subio_unfiltered_wo_names.txt", header = TRUE, sep = "\t")
cts <- cts[rowSums(cts) != 0,]
cts <- as.data.frame(t(cts))
cts <- rownames_to_column(cts, var = "V1")
cts$V1 <- c("PTC", "PTC", "PC", "PC", "aCD38", "aCD38")
cts$V1 <- factor(cts$V1, levels = c("PC", "PTC", "aCD38"))

mtx <- read.delim(file = "data/Bulk_fibro_subio_unfiltered_wo_names.txt", header = TRUE, sep = "\t")
mtx <- mtx[rowSums(mtx) != 0,]
pca <- prcomp(t(mtx), scale. = TRUE)
autoplot(pca, data = cts, colour = "V1", size = 3) +
  theme_bw()

