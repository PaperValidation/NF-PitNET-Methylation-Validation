#!/usr/bin/env Rscript
# install missing lightweight R packages (Bioconductor + CRAN)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("survminer", "timeROC", "rms"), ask = FALSE, update = FALSE)
# ChAMP is large; install only if network allows
# BiocManager::install("ChAMP", ask=FALSE, update=FALSE)
cat("=== post-install check ===\n")
for (p in c("survminer","timeROC","rms","ChAMP")) {
  cat(sprintf("%-12s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "OK", "MISSING")))
}
