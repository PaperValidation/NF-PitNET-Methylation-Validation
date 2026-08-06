#!/usr/bin/env Rscript
# ============================================================================
# Møller 2025 (PMID 41369112) raw EPIC 850K .idat -> beta matrix (cg probes)
# ----------------------------------------------------------------------------
# Input  (user-supplied, from Mendeley Data DOI 10.17632/vr2wb38r6t.1, on DATASET_PATH G: drive):
#   G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Mendeley_Moller2025/Deidentified_IDATs/   # .idat
#   G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Mendeley_Moller2025/SampleSheet.xlsx      # Sentrix_ID, Sentrix_Position, Sample_ID(FFPE-MM-X)
# Output (written alongside the raw data on G:, NOT the C: project dir):
#   <Mendeley_Moller2025>/mendeley_beta_cg.rds     # matrix: cg probes x samples, colnames = IDAT basename
#   <Mendeley_Moller2025>/mendeley_beta_cg.csv.gz  # same, gzipped TSV (colnames = IDAT basename)
#   <Mendeley_Moller2025>/mendeley_sample_manifest.csv  # basename(IDAT) -> Sample_ID(FFPE-MM-X), Sentrix_ID, Sentrix_Position
#
# NOTE on sample identity / linkage key:
#   The Mendeley SampleSheet uses Sample_ID = "FFPE-MM-X", but Møller's
#   Supplementary_File_S1.xlsx uses an INTEGER Sample_ID and an "IDAT" column
#   equal to "<Sentrix_ID>_<Sentrix_Position>". The idat FILE basename is also
#   "<Sentrix_ID>_<Sentrix_Position>". So we KEEP the idat basename as the
#   matrix column name (it equals S1's IDAT key) and emit a manifest that
#   maps basename <-> FFPE-MM-X for reporting. The strict-projection script
#   links to S1 via the IDAT basename, not via the FFPE-MM-X string.
#
# Run in WSL (matches the Neou idat pipeline):
#   WSL$ cd /mnt/c/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET
#   WSL$ Rscript scripts/moller_idat_to_beta.R
#
# Requires: BiocManager::install(c("minfi","readxl","IlluminaHumanMethylationEPICmanifest"))
# Memory: 117 EPIC 850K arrays -> RGset a few GB RAM; run on >= 16 GB.
# ============================================================================
suppressMessages({
  library(minfi)
  library(readxl)
  library(dplyr)
})

# Raw Mendeley data lives on the DATASET_PATH G: drive (NOT the C: project dir)
moller_root_wsl <- "/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Mendeley_Moller2025"
moller_root_win <- "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Mendeley_Moller2025"
root     <- if (.Platform$OS.type == "windows") moller_root_win else moller_root_wsl
idat_dir <- file.path(root, "Deidentified_IDATs")
sheet    <- file.path(root, "SampleSheet.xlsx")
out_rds  <- file.path(root, "mendeley_beta_cg.rds")
out_csv  <- file.path(root, "mendeley_beta_cg.csv.gz")
manifest <- file.path(root, "mendeley_sample_manifest.csv")

stopifnot(dir.exists(idat_dir), file.exists(sheet))

# ---- read idat (EPIC 850K; minor per-file probe-count variance -> force=TRUE) ----
# idat files live in Deidentified_IDATs/<Sentrix_ID>/<Sentrix_ID>_<Pos>[_Grn/_Red].idat
ss <- readxl::read_excel(sheet)
ss <- ss %>% dplyr::mutate(basename = paste(.data$Sentrix_ID, .data$Sentrix_Position, sep = "_"))
bns <- file.path(idat_dir, ss$Sentrix_ID, paste0(ss$Sentrix_ID, "_", ss$Sentrix_Position))
miss <- sum(!file.exists(paste0(bns, "_Grn.idat")))
if (miss > 0) warning(sprintf("%d SampleSheet arrays missing Grn idat on disk", miss))

cat("[1] Reading", length(bns), "EPIC 850K arrays (force=TRUE) ...\n")
rg   <- minfi::read.metharray(basenames = bns, force = TRUE)
beta <- minfi::getBeta(rg)
colnames(beta) <- basename(bns)   # short basename = S1 IDAT key
cat("    raw beta:", nrow(beta), "probes x", ncol(beta), "arrays\n")

# keep only cg (CpG) probes; drop rs/snp/ch probes
cg <- beta[grepl("^cg", rownames(beta)), , drop = FALSE]
cat("    cg probes retained:", nrow(cg), "\n")

# ---- linkage: idat basename = "<Sentrix_ID>_<Sentrix_Position>" ----
cat("[2] SampleSheet cols:", paste(colnames(ss), collapse = ", "), "\n")
lut <- setNames(as.character(ss$Sample_ID), ss$basename)   # basename -> FFPE-MM-X

present <- colnames(cg)
keep    <- present[present %in% names(lut)]
cg      <- cg[, keep, drop = FALSE]
cat("    matched arrays kept:", ncol(cg), "of", length(lut), "SampleSheet entries\n")
if (length(keep) == 0) stop("No idat basenames matched the SampleSheet — check Sentrix_ID/Position.")

# colnames stay as the IDAT basename (= S1 IDAT key, universal join key)
# emit manifest: basename(IDAT) <-> Sample_ID(FFPE-MM-X) + Sentrix info
man <- ss %>% dplyr::select(basename, Sample_ID, Sentrix_ID, Sentrix_Position) %>%
  dplyr::filter(basename %in% keep) %>%
  dplyr::arrange(basename)
write.csv(man, manifest, row.names = FALSE)
cat("    wrote manifest ->", manifest, "(", nrow(man), "rows )\n")

# ---- save ----
saveRDS(cg, out_rds)
cat("[3] Saved rds ->", out_rds, "(", nrow(cg), "x", ncol(cg), ")\n")

gz <- gzfile(out_csv, "wt")
write.table(cg, file = gz, sep = "\t", quote = FALSE, col.names = NA)
close(gz)
cat("    Saved csv.gz ->", out_csv, "\n")
cat("DONE.\n")
