#!/usr/bin/env Rscript
# Export Neou EPIC beta matrix (cg-keyed) from RDS to a Python-readable binary.
# Avoids rpy2: writes dims + float64 matrix + probe/sample name files.
args <- commandArgs(trailingOnly = TRUE)
RDS <- "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/E-MTAB-7762/beta_out/neou_beta_cg.rds"
OUTDIR <- "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/E-MTAB-7762/beta_out"
if (length(args) >= 1) RDS <- args[1]
if (length(args) >= 2) OUTDIR <- args[2]

cat("[1] Loading", RDS, "\n")
B <- readRDS(RDS)
cat("    dims:", nrow(B), "probes x", ncol(B), "samples\n")
cat("    rownames head:", head(rownames(B), 3), "\n")
stopifnot(grepl("^cg", rownames(B)[1]))   # must be cg-keyed

bin <- file.path(OUTDIR, "neou_beta_full.bin")
cat("[2] Writing binary ->", bin, "\n")
con <- file(bin, "wb")
writeBin(as.integer(nrow(B)), con, size = 4)
writeBin(as.integer(ncol(B)), con, size = 4)
writeBin(as.numeric(B), con, size = 8)    # float64
close(con)
cat("    bytes:", file.info(bin)$size, "\n")

cat("[3] Writing probe + sample name files\n")
writeLines(rownames(B), file.path(OUTDIR, "neou_beta_probes.txt"))
writeLines(colnames(B), file.path(OUTDIR, "neou_beta_samples.txt"))
cat("[done]\n")
