#!/usr/bin/env Rscript
# ============================================================================
# Normalization-sensitivity analysis (Neou E-MTAB-7762, EPIC 850K)
# ----------------------------------------------------------------------------
# Reprocesses the Neou idat with four normalizations (preprocessRaw / Noob /
# Funnorm / Quantile) and re-runs the nearest-centroid projection to the
# GSE283928 c0/c2/c3/c4 centroids, quantifying how stable the c0 assignment is.
# Feeds manuscript Table 3 / Section 3.1.
#
# Provenance: recovered from the archived WSL working copy (E:\Data Temp\
# NF-PitNET\_wsl_scripts\norm_sensitivity.R, 2026-08-09) and adapted to the
# shared path loader (scripts/paths.R) + DATASET_PATH_TEMP output convention.
#
# Inputs (raw, DATASET_PATH G:):
#   E-MTAB-7762/idat_flat/                    (Neou .idat pairs)
#   Moller2025_derived/centroid_c0_c2_c3_c4.csv
# Outputs (intermediate, DATASET_PATH_TEMP E:):
#   intermediates/norm_sensitivity_assignment.rds
# Run under WSL: Rscript scripts/norm_sensitivity.R
# ============================================================================
suppressMessages({library(minfi)})
source("C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/scripts/paths.R")
BASE <- paths$DATASET_PATH_WSL          # /mnt/e/... (raw data)
TEMP <- paths$DATASET_PATH_TEMP_WSL      # /mnt/e/Data Temp/NF-PitNET (intermediates)

idat.dir <- file.path(BASE, "E-MTAB-7762", "idat_flat")
grn <- list.files(idat.dir, pattern="_Grn.idat$", full.names=FALSE)
basenames <- sub("_Grn.idat$","", grn)
message("found ", length(basenames), " idat pairs")
rg <- read.metharray(file.path(idat.dir, basenames), force=TRUE, verbose=FALSE)
message("RGChannelSet dims: ", nrow(rg), " x ", ncol(rg))

beta_raw   <- getBeta(preprocessRaw(rg))
beta_noob  <- tryCatch(getBeta(preprocessNoob(rg)), error=function(e){message("noob failed: ",e$message); NULL})
beta_fun   <- tryCatch(getBeta(preprocessFunnorm(rg, verbose=FALSE)), error=function(e){message("funnorm failed: ",e$message); NULL})
beta_quant <- tryCatch(getBeta(preprocessQuantile(rg)), error=function(e){message("quantile failed: ",e$message); NULL})

cen <- read.csv(file.path(BASE, "Moller2025_derived", "centroid_c0_c2_c3_c4.csv"))
cenprobes <- as.character(cen$probe)
common <- intersect(rownames(beta_raw), cenprobes)
message("common probes with centroid: ", length(common))
C <- as.matrix(cen[match(common, cen$probe), c("c0","c2","c3","c4")]); rownames(C) <- common

assign_to <- function(B){
  Bc <- B[common, , drop=FALSE]
  d0 <- colSums((Bc - C[,"c0"])^2); d2 <- colSums((Bc - C[,"c2"])^2)
  d3 <- colSums((Bc - C[,"c3"])^2); d4 <- colSums((Bc - C[,"c4"])^2)
  apply(rbind(d0,d2,d3,d4), 2, which.min)
}
res <- list()
for (nm in c("raw","noob","fun","quant")) {
  B <- get(paste0("beta_", nm)); if (is.null(B)) next
  res[[nm]] <- assign_to(B)
  message(nm, " c0 fraction: ", round(mean(res[[nm]]==1),3))
}
dir.create(file.path(TEMP, "intermediates"), showWarnings = FALSE, recursive = TRUE)
saveRDS(res, file.path(TEMP, "intermediates", "norm_sensitivity_assignment.rds"))
if (!is.null(res$raw)) {
  ag <- sapply(names(res), function(nm) mean(res[[nm]]==res$raw))
  cat("agreement with preprocessRaw assignment:\n"); print(round(ag,4))
}
message("DONE")
