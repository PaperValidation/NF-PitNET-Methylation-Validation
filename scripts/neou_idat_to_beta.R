#!/usr/bin/env Rscript
# ============================================================================
# Neou 2020 (E-MTAB-7762) raw .idat -> beta matrix (Illumina 450K)
# Run under WSL: Rscript neou_idat_to_beta.R
# Reads idat from /mnt/g .../E-MTAB-7762, outputs beta + sample annotation
# ============================================================================

suppressMessages({
  library(minfi)
  library(data.table)
})

BASE <- "/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/E-MTAB-7762"
SDRF <- file.path(BASE, "E-MTAB-7762.sdrf.txt")
OUT  <- file.path(BASE, "beta_out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("== base dir exists:", dir.exists(BASE), "\n")
cat("== listing idat files ==\n")
idat <- list.files(BASE, pattern = "\\.idat(\\.gz)?$", recursive = TRUE, full.names = TRUE)
cat(sprintf("found %d idat files\n", length(idat)))
if (length(idat) == 0) stop("no idat found -- unzip the raw.*.zip first")

# ---- 1. read RGset -------------------------------------------------------
cat("== read.metharray (RGset) ==\n")
t0 <- Sys.time()
idat_files <- unlist(lapply(paste0(BASE, "/raw_", 1:5),
                              function(d) list.files(d, pattern = "\\.idat(\\.gz)?$",
                                                    full.names = TRUE)))
cat(sprintf("enumerated %d idat files (raw_1..5 only)\n", length(idat_files)))
# The 5 raw.*.zip split each sample's Grn/Red into DIFFERENT subdirs, so a
# basename-with-dir cannot locate both channels. Co-locate every idat into a
# single flat symlink dir, then pass the flat sample stems (vector form) to
# read.metharray, which appends _Grn.idat/_Red.idat itself.
flat <- file.path(BASE, "idat_flat")
if (dir.exists(flat)) unlink(flat, recursive = TRUE)
dir.create(flat, recursive = TRUE)
for (f in idat_files) file.symlink(f, file.path(flat, basename(f)))
flat_idat <- list.files(flat, pattern = "\\.idat(\\.gz)?$", full.names = TRUE)
stems <- sub("_Red\\.idat.*$", "", sub("_Grn\\.idat.*$", "", flat_idat))
stems <- unique(stems)
cat(sprintf("flat co-located idat: %d | unique sample stems: %d\n",
            length(flat_idat), length(stems)))
RGset <- read.metharray(stems, verbose = TRUE, force = TRUE)
cat(sprintf("RGset read in %.1f min; samples=%d probes=%d\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            ncol(RGset), nrow(RGset)))
cat("== RGset rownames sample ==\n")
print(head(rownames(RGset), 6))
cat(sprintf("nProbes RGset: %d ; cg-like rownames: %d\n",
            nrow(RGset), sum(grepl("^cg", rownames(RGset)))))

# ---- 2. preprocess + beta ------------------------------------------------
cat("== preprocessRaw + getBeta ==\n")
MSet <- preprocessRaw(RGset)
beta <- getBeta(MSet)
cat(sprintf("beta dims: %d probes x %d samples\n", nrow(beta), ncol(beta)))

# ---- 3. detection QC (flag poor samples) --------------------------------
cat("== detectionP QC ==\n")
detP <- detectionP(RGset)
qc <- data.frame(
  sample = colnames(beta),
  n_fail_detP_gt_0.01 = colSums(detP > 0.01, na.rm = TRUE),
  mean_beta = round(colMeans(beta, na.rm = TRUE), 4),
  stringsAsFactors = FALSE
)
print(head(qc[order(-qc$n_fail_detP_gt_0.01), ], 10))

# ---- 4. map numeric rownames -> cg + annotate (EPIC hg19) --------------
cat("== map probe IDs to cg names (EPIC) ==\n")
anno <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
mapA <- setNames(rownames(anno), anno$AddressA)
mapB <- setNames(rownames(anno), anno$AddressB)
rid <- rownames(beta)
cg <- ifelse(!is.na(mapA[rid]), mapA[rid], mapB[rid])
mapped <- !is.na(cg)
cat(sprintf("beta probes=%d ; mapped to cg=%d ; unmapped=%d\n",
            nrow(beta), sum(mapped), sum(!mapped)))
beta_cg <- beta[mapped, , drop = FALSE]
rownames(beta_cg) <- cg[mapped]
common <- intersect(rownames(beta_cg), rownames(anno))
beta_a <- beta_cg[common, , drop = FALSE]
anno_a <- as.data.frame(anno[common, c("chr", "pos", "UCSC_RefGene_Name",
                                       "UCSC_RefGene_Group", "Relation_to_Island")])
anno_a <- anno_a[rownames(beta_a), ]
rownames(anno_a) <- rownames(beta_a)
cat(sprintf("annotated beta: %d probes (cg-keyed)\n", nrow(beta_a)))

# ---- 5. sample annotation from SDRF -------------------------------------
cat("== parse SDRF for sample annotation ==\n")
sd <- fread(SDRF, fill = TRUE, sep = "\t", header = TRUE, check.names = FALSE)
adf_col <- colnames(sd)[grepl("array data file", tolower(colnames(sd)))][1]
ct_col  <- colnames(sd)[grepl("characteristics\\[cell type\\]", tolower(colnames(sd)))][1]
ds_col  <- colnames(sd)[grepl("characteristics\\[disease staging\\]", tolower(colnames(sd)))][1]
ind_col <- colnames(sd)[grepl("characteristics\\[individual\\]", tolower(colnames(sd)))][1]

sd[[".adf"]]  <- basename(as.character(sd[[adf_col]]))
sd[[".stem"]] <- sub("_Grn.*|_Red.*", "", sd[[".adf"]])
ann <- unique(sd[, c(ind_col, ct_col, ds_col, ".stem"), with = FALSE])
setnames(ann, c("individual", "cell_type", "disease_staging", "stem"))
# NFPT = non-secreting (Silent or Gonadotrophin staging)
ann[["nfpt"]] <- ifelse(ann[["disease_staging"]] %in% c("Silent", "Gonadotrophin"), 1L, 0L)
ann[["secretion"]] <- ifelse(ann[["nfpt"]] == 1L, "non-secreting", "secreting")
# lineage mapping (pituitary transcription-factor lineage)
lin <- c("gonadtroph"="SF1", "Null-cell"="SF1",
         "corticotroph"="TPIT",
         "somatotroph"="PIT1", "mammotroph"="PIT1", "thyrotroph"="PIT1",
         "Plurihormonal_PIT-1"="PIT1", "Mixed_GH-PRL"="PIT1")
ann[["lineage"]] <- lin[ann[["cell_type"]]]
# align annotation rows to beta column order
ann <- ann[match(colnames(beta), ann[["stem"]]), ]
rownames(ann) <- ann[["stem"]]
cat(sprintf("annotation rows=%d  matches beta cols=%s\n",
            nrow(ann), nrow(ann) == ncol(beta)))
cat("NFPT count:", sum(ann[["nfpt"]]), " | lineage dist:\n")
print(table(ann[["lineage"]], useNA = "ifany"))

# ---- 6. save -------------------------------------------------------------
cat("== saving outputs ==\n")
saveRDS(beta,    file.path(OUT, "neou_beta_raw.rds"))
saveRDS(beta_cg, file.path(OUT, "neou_beta_cg.rds"))
saveRDS(beta_a,  file.path(OUT, "neou_beta_annotated.rds"))
fwrite(as.data.frame(beta_a), file.path(OUT, "neou_beta_annotated.tsv.gz"),
       sep = "\t", compress = "gzip", row.names = TRUE, quote = FALSE)
fwrite(as.data.frame(anno_a), file.path(OUT, "neou_probe_annotation.tsv.gz"),
       sep = "\t", compress = "gzip", row.names = TRUE, quote = FALSE)
fwrite(qc, file.path(OUT, "neou_sample_qc.csv"), sep = ",")
fwrite(ann, file.path(OUT, "neou_sample_annotation.csv"), sep = ",")
cat("== DONE. outputs in", OUT, "\n")
cat("  neou_beta_annotated.tsv.gz :", nrow(beta_a), "probes x", ncol(beta_a), "samples\n")
cat("  neou_sample_annotation.csv :", nrow(ann), "samples (", sum(ann[["nfpt"]]), "NFPT )\n")
