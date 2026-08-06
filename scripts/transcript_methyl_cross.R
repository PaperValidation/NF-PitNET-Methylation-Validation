suppressMessages({
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, quietly = TRUE)
  library(data.table, quietly = TRUE)
})

BASE <- "/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
RES  <- file.path(BASE, "results", "GSE209903")
OUT  <- file.path(BASE, "results", "risk_model")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## GSE209903 NFPT-up genes (NFPT higher = log2FC < 0)
de <- fread(file.path(RES, "FPT_vs_NFPT_full.csv"))
nfpt_up <- de[padj < 0.05 & log2FoldChange < -1 & baseMean > 50 & !is.na(SYMBOL) & SYMBOL != "", SYMBOL]
cat("NFPT-up genes (GSE209903):", length(nfpt_up), "\n")

## Map to EPIC probes (gene -> probes)
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
gcol <- grep("UCSC_RefGene_Name", colnames(anno), value = TRUE)[1]
anno_dt <- data.table(Probe = rownames(anno), Gene = anno[[gcol]])
nfpt_up_set <- unique(nfpt_up)
mapped <- anno_dt[Gene %in% nfpt_up_set]
cat("NFPT-up genes mapped to EPIC probes:", length(unique(mapped$Gene)),
    "/", length(nfpt_up_set), "\n")
cat("Total EPIC probes for NFPT-up program:", nrow(mapped), "\n")

## Save probe list for Python-side per-cluster methylation test in GSE283928
fwrite(mapped[, .(Probe, Gene)], file.path(OUT, "nfpt_up_epic_probes.csv"))
cat("Saved nfpt_up_epic_probes.csv\n")

## Also build a generic "NFPT transcriptional program" summary
summary <- data.frame(
  n_nfpt_up_genes = length(nfpt_up_set),
  n_mapped_to_epic = length(unique(mapped$Gene)),
  n_epic_probes = nrow(mapped)
)
fwrite(summary, file.path(OUT, "crossmap_summary.csv"))
print(summary)
