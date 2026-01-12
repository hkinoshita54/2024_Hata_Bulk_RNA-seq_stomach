# result table
description <- paste0(test_group, "_vs_", ref_group)
res <- results(dds, contrast = c("exp_group", test_group, ref_group)) %>% data.frame()
res <- rownames_to_column(res, var = "ensmusg")
res <- left_join(res, anno_mm, by = c("ensmusg" = "ensembl_gene_id"))
openxlsx2::write_xlsx(res, file.path(res_path, paste0(description, ".xlsx")))

## volcano plot: WIP ----
# library(EnhancedVolcano)
# keyvals.colour <- ifelse(res$padj > 0.01 | abs(res$log2FoldChange) < 2, 'gray', 'orange')
# names(keyvals.colour)[keyvals.colour == 'gray'] <- 'other'
# names(keyvals.colour)[keyvals.colour == 'orange'] <- '|logFC|>2 & padj<0.01'
# EnhancedVolcano(res, lab = res$external_gene_name, x = 'log2FoldChange', y = 'pvalue',
#                 # xlim = c(-10, 10), ylim = c(0, -log10(10e-100)), legendPosition = "right",
#                 # selectLab = c('Muc5ac','Gkn2'),
#                 drawConnectors = TRUE, widthConnectors = 0.5, arrowheads = FALSE,
#                 # pCutoff = 10e-6, FCcutoff = 2, col = c('gray', 'gray', 'gray', 'orange'),
#                 title = paste0(test_group, " vs ", ref_group), subtitle = NULL, caption = NULL,
#                 colCustom = keyvals.colour,
#                 pointSize = 1.0, labSize = 4.0, cutoffLineWidth = 0.0, colAlpha = 1)
# ggsave(file.path("plot", paste0("volcano_", description, ".png")), height = 6, width = 5, dpi = 150)
# ----

# fgsea

## make gene rank
rank <- res %>% dplyr::select(external_gene_name, stat) %>% na.omit() %>% 
  group_by(external_gene_name) %>% summarize(stat=mean(stat)) %>% deframe()

## HALLMARK with plots
fgseaRes <- fgsea(pathways = collections$H, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_H_", description, ".xlsx")))

fgseaRes$pathway <- sub("HALLMARK_", "", fgseaRes$pathway)
fgseaResUp <- fgseaRes[fgseaRes$padj<0.25 & NES>0,]
fgseaResDn <- fgseaRes[fgseaRes$padj<0.25 & NES<0,]
cols <- brewer.pal(3, "Set1")

ggplot(fgseaResUp, aes(x = NES, y = reorder(pathway, NES))) +
  geom_col(fill = cols[1]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  # scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(test_group, " vs ", ref_group, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Up_", description, ".png"), path = plot_path, width = 3.5, height = 0.5 + nrow(fgseaResUp)/8, units = "in", dpi = 300) 

ggplot(fgseaResDn, aes(x = NES, y = reorder(pathway, NES, decreasing = TRUE))) +
  geom_col(fill = cols[2]) +
  scale_x_continuous(expand=c(0,0), position = "top") +
  scale_y_discrete(position = "right") +
  labs(x="NES", y=NULL, title=paste0(test_group, " vs ", ref_group, " (FDR<0.25)")) +
  theme_classic() +
  theme(plot.title = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"))
ggsave(paste0("GSEA_H_Dn_", description, ".png"), path = plot_path, width = 3.5, height = 0.5 + nrow(fgseaResDn)/8, units = "in", dpi = 300)

# KEGG, REACTOME, BP, C6 with collapsePathways (without plotting)
run_fgsea <- function(gene_set, gene_set_name){
  fgseaRes <- fgsea(pathways = gene_set, stats = rank, eps=0.0, minSize = 10, maxSize = 500)
  fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.25], gene_set, rank)
  fgseaRes <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
  openxlsx2::write_xlsx(fgseaRes, file.path(res_path, paste0("GSEA_", gene_set_name, "_", description, ".xlsx")))
}

for (i in 2:5){
  run_fgsea(collections[[i]], names(collections)[i])
}