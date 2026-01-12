library(GCclassifier)
data("GSE62254")
data("GSE62254_subtype")
GSE62254[1:5, 1:5]
identical(colnames(GSE62254), GSE62254_subtype$GEO_ID)

emp.res <- classifyGC(
  Expr = GSE62254, ## gene expression profile with log2 transformation
  method = 'EMP', ## subtyping system
  idType = 'SYMBOL', ## the gene identifier type in gene expression profile
  useMinPosterior = F ## whether the minPosterior threshold is also used for EMT subtyping
)
table(emp.res$subtype, GSE62254_subtype$Subgroup)
