#!/usr/bin/env Rscript
# Finalize Neou beta: annotate (cg-keyed) + sample annotation, from saved RDS.
# Avoids re-reading idat (5 min). Run under WSL.
suppressMessages({ library(minfi); library(data.table) })
BASE   <- "/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/E-MTAB-7762"
OUT    <- file.path(BASE, "beta_out")
SDRF   <- file.path(BASE, "E-MTAB-7762.sdrf.txt")

cat("== loading saved beta_raw.rds ==\n")
beta <- readRDS(file.path(OUT, "neou_beta_raw.rds"))
cat(sprintf("beta dims: %d probes x %d samples\n", nrow(beta), ncol(beta)))
rid <- rownames(beta)
cat("rownames head:", head(rid, 3), "| cg-like:", sum(grepl("^cg", rid)), "\n")

cat("== annotate (EPIC hg19) ==\n")
anno <- getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
if (!grepl("^cg", rid[1])) {
  mapA <- setNames(rownames(anno), anno$AddressA)
  mapB <- setNames(rownames(anno), anno$AddressB)
  cg <- ifelse(!is.na(mapA[rid]), mapA[rid], mapB[rid])
  mapped <- !is.na(cg)
  beta <- beta[mapped, , drop = FALSE]; rownames(beta) <- cg[mapped]; rid <- rownames(beta)
  cat(sprintf("mapped numeric->cg: %d probes\n", sum(mapped)))
}
common <- intersect(rid, rownames(anno))
cat(sprintf("cg probes overlapping EPIC anno: %d / %d\n", length(common), nrow(beta)))
beta_a <- beta[common, , drop = FALSE]
anno_a <- as.data.frame(anno[common, c("chr", "pos", "UCSC_RefGene_Name",
                                       "UCSC_RefGene_Group", "Relation_to_Island")])
anno_a <- anno_a[rownames(beta_a), ]; rownames(anno_a) <- rownames(beta_a)

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
cat("== disease_staging raw values (in ann) ==\n")
print(table(ann[["disease_staging"]]))
ann[["nfpt"]] <- ifelse(ann[["disease_staging"]] %in% c("Silent", "Gonadotrophin"), 1L, 0L)
ann[["secretion"]] <- ifelse(ann[["nfpt"]] == 1L, "non-secreting", "secreting")
lin <- c("gonadtroph"="SF1", "Null-cell"="SF1",
         "corticotroph"="TPIT",
         "somatotroph"="PIT1", "mammotroph"="PIT1", "thyrotroph"="PIT1",
         "Plurihormonal_PIT-1"="PIT1", "Mixed_GH-PRL"="PIT1")
ann[["lineage"]] <- lin[ann[["cell_type"]]]
ann <- ann[match(colnames(beta), ann[["stem"]]), ]; rownames(ann) <- ann[["stem"]]
cat(sprintf("annotation rows=%d matches beta cols=%s ; NFPT=%d\n",
            nrow(ann), nrow(ann) == ncol(beta), sum(ann[["nfpt"]])))

cat("== saving outputs ==\n")
saveRDS(beta,    file.path(OUT, "neou_beta_cg.rds"))
saveRDS(beta_a,  file.path(OUT, "neou_beta_annotated.rds"))
fwrite(as.data.frame(beta_a), file.path(OUT, "neou_beta_annotated.tsv.gz"),
       sep = "\t", compress = "gzip", row.names = TRUE, quote = FALSE)
fwrite(as.data.frame(anno_a), file.path(OUT, "neou_probe_annotation.tsv.gz"),
       sep = "\t", compress = "gzip", row.names = TRUE, quote = FALSE)
fwrite(ann, file.path(OUT, "neou_sample_annotation.csv"), sep = ",")
cat("== DONE ==\n")
cat("  neou_beta_cg.rds           :", nrow(beta), "probes x", ncol(beta), "samples\n")
cat("  neou_beta_annotated.tsv.gz:", nrow(beta_a), "probes x", ncol(beta_a), "samples\n")
cat("  neou_sample_annotation.csv :", nrow(ann), "samples (", sum(ann[["nfpt"]]), "NFPT )\n")
