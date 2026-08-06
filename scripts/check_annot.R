#!/usr/bin/env Rscript
suppressMessages(library(data.table))
m <- fread("/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/GSE51618_series_matrix.txt.gz",
           sep="\t", header=TRUE, skip=59, fill=TRUE, data.table=FALSE)
cat("First 6 colnames:\n")
print(colnames(m)[1:6])
cat("\nHead of first 3 cols:\n")
print(head(m[,1:3]))
# check if hgu133a.db available
if (requireNamespace("hgu133a.db", quietly=TRUE)) {
  cat("\nhgu133a.db: AVAILABLE\n")
} else {
  cat("\nhgu133a.db: MISSING\n")
}
