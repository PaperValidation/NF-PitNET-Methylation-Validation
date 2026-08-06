#!/usr/bin/env Rscript
# check installed R packages relevant to the NF-PitNET pipeline
pkgs <- c("GEOquery","limma","edgeR","DESeq2","sva","WGCNA","glmnet",
          "survival","survminer","clusterProfiler","fgsea","Seurat",
          "SingleR","minfi","ChAMP","DMRcate","tidyverse","data.table",
          "ggplot2","pROC","timeROC","rms","AnnotationHub","ExperimentHub")
cat("R version:", as.character(getRversion()), "\n")
inst <- installed.packages()[,1]
for (p in pkgs) {
  cat(sprintf("%-18s %s\n", p, ifelse(p %in% inst, "OK", "MISSING")))
}
