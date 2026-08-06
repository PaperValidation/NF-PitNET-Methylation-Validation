#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# install lightweight pure-R packages only (avoid Deriv/nloptr chain from rms)
for (p in c("survminer", "timeROC", "cmprsk")) {
  cat("=== installing", p, "===\n")
  BiocManager::install(p, ask = FALSE, update = FALSE)
}
cat("=== check ===\n")
for (p in c("survminer","timeROC","cmprsk","rms")) {
  cat(sprintf("%-12s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "OK", "MISSING")))
}
