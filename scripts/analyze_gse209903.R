#!/usr/bin/env Rscript
# GSE209903: FPT vs NFPT precomputed DESeq2 results (Silva-Junior 2023 paired RNA-seq)
# Extract significant DE genes, annotate ENSG->symbol, cross-check vs Moller 5-cluster key genes.
suppressMessages({ library(data.table); library(dplyr); library(ggplot2) })
BASE <- "/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
DE   <- file.path(BASE, "GSE209903_FPT_vs_NFPT.txt.gz")
OUT  <- file.path(BASE, "results", "GSE209903")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dt <- fread(DE, sep = "\t", header = TRUE, data.table = FALSE)
cat("DE rows:", nrow(dt), "\n")
dt$sig <- dt$padj < 0.05 & abs(dt$log2FoldChange) > 1
cat("Significant (padj<0.05 & |log2FC|>1):", sum(dt$sig), "\n")

## ENSG -> symbol via org.Hs.eg.db if available, else keep ENSG
sym <- dt$ID
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressMessages({
    library(org.Hs.eg.db)
    map <- mapIds(org.Hs.eg.db, keys = dt$ID, column = "SYMBOL", keytype = "ENSEMBL")
  })
  sym <- ifelse(is.na(map), dt$ID, map)
}
dt$SYMBOL <- sym

## top tables
up   <- dt %>% filter(sig, log2FoldChange > 0) %>% arrange(padj)
down <- dt %>% filter(sig, log2FoldChange < 0) %>% arrange(padj)
fwrite(dt, file.path(OUT, "FPT_vs_NFPT_full.csv"))
fwrite(head(up, 200),   file.path(OUT, "NFPT_up_top200.csv"))
fwrite(head(down, 200), file.path(OUT, "NFPT_down_top200.csv"))

cat("\n=== TOP 20 UP in NFPT (FPT<NFPT means negative log2FC => gene higher in NFPT) ===\n")
# Note: FPT_vs_NFPT => positive log2FC = higher in FPT. So NFPT-up = negative log2FC.
nfpt_up <- dt %>% filter(sig, log2FoldChange < 0) %>% arrange(log2FoldChange)
print(head(nfpt_up[, c("SYMBOL","log2FoldChange","padj")], 20))
cat("\n=== TOP 20 UP in FPT ===\n")
fpt_up <- dt %>% filter(sig, log2FoldChange > 0) %>% arrange(desc(log2FoldChange))
print(head(fpt_up[, c("SYMBOL","log2FoldChange","padj")], 20))

## cross-check with Moller 2025 key genes (5-cluster signature, from full text)
moller_key <- c("CCDC141","ATP1A3","ANK2","CASC15","HDC","NCKAP5","GUCA1C","ATP2C1",
                "SF1","NR5A1","GATA3","TPIT","PIT1","TBX19","PROP1","POU1F1")
cat("\n=== Moller key genes present in this DE set ===\n")
mm <- dt %>% filter(SYMBOL %in% moller_key)
print(mm[, c("SYMBOL","log2FoldChange","padj","sig")])

## volcano
png(file.path(OUT, "volcano.png"), width = 900, height = 700)
with(dt, {
  col <- ifelse(sig & log2FoldChange > 0, "red",
         ifelse(sig & log2FoldChange < 0, "blue", "grey"))
  plot(log2FoldChange, -log10(padj), col = col, pch = 19, cex = 0.5,
       xlab = "log2FC (FPT / NFPT)", ylab = "-log10 padj)",
       main = "GSE209903: FPT vs NFPT (65 RNA-seq)")
  abline(h = -log10(0.05), v = c(-1,1), lty = 2, col = "darkgrey")
})
dev.off()
cat("\nSaved to", OUT, "\n")
