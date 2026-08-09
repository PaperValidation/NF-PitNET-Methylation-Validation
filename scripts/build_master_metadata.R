#!/usr/bin/env Rscript
source("C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/scripts/paths.R")
# Build a unified master sample metadata table across downloaded GEO cohorts,
# filtering to non-functioning pituitary neuroendocrine tumors (NFPT/NFPA)
# where applicable, and standardizing phenotype columns.
suppressMessages({
  library(dplyr); library(stringr); library(readr)
})
BASE <- paths$DATASET_PATH_WSL
SAMP <- file.path(paths$DATASET_PATH_TEMP_WSL, "samples")
OUT  <- file.path(paths$DATASET_PATH_TEMP_WSL, "metadata")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- GSE51618: explicit invasion labels, microarray ----
g1 <- read_csv(file.path(SAMP, "GSE51618_samples.csv"), show_col_types = FALSE) %>%
  mutate(
    cohort = "GSE51618",
    geo_sample = sample_geo_accession,
    platform = sample_platform_id,
    title = sample_title,
    source = sample_source_name_ch1,
    characteristics = sample_characteristics_ch1,
    # classify by title keywords
    group = case_when(
      str_detect(title, "non-invasion") ~ "NFPA_non-invasive",
      str_detect(title, "invasion")     ~ "NFPA_invasive",
      str_detect(title, "pituitary gland") ~ "normal_pituitary",
      TRUE ~ "other"
    ),
    nfpt = group %in% c("NFPA_non-invasive","NFPA_invasive"),
    omics = "bulk_transcriptome_microarray"
  ) %>%
  select(cohort, geo_sample, platform, title, source, characteristics, group, nfpt, omics)

## ---- GSE207937: PitNET spectrum (EPIC methylation), subtype annotated ----
## Silva-Junior 2023 subtypes are GH/ACTH/NS/NA only; NFPT not explicitly tagged.
## We keep subtype verbatim; NS = non-secretory (NFPT candidate), NA = unknown.
g2 <- read_csv(file.path(SAMP, "GSE207937_samples.csv"), show_col_types = FALSE) %>%
  mutate(
    cohort = "GSE207937",
    geo_sample = sample_geo_accession,
    platform = sample_platform_id,
    title = sample_title,
    source = sample_source_name_ch1,
    characteristics = sample_characteristics_ch1,
    subtype = str_trim(str_extract(characteristics, "(?<=subtypes:\\s)\\S+")),
    nfpt_candidate = subtype == "NS",         # non-secretory ~ NFPT
    group = case_when(
      subtype == "GH"   ~ "GH_PitNET",
      subtype == "ACTH" ~ "ACTH_PitNET",
      subtype == "NS"   ~ "NFPT_candidate",
      TRUE              ~ "unknown_subtype"
    ),
    nfpt = NA,                                  # not reliably derivable from GEO meta
    omics = "bulk_methylation_EPIC"
  ) %>%
  select(cohort, geo_sample, platform, title, source, characteristics, subtype, group, nfpt_candidate, nfpt, omics)

## ---- GSE209903: matched RNA-seq to GSE207937, filter NFPT by transcriptome_group ----
g3 <- read_csv(file.path(SAMP, "GSE209903_samples.csv"), show_col_types = FALSE) %>%
  mutate(
    cohort = "GSE209903",
    geo_sample = sample_geo_accession,
    platform = sample_platform_id,
    title = sample_title,
    source = sample_source_name_ch1,
    characteristics = sample_characteristics_ch1,
    tgroup = str_extract(characteristics, "(?<=transcriptome_group:\\s)\\S+"),
    # functional vs non-functional pituitary tumor
    nfpt = str_detect(characteristics, regex("non-function|FPT_nf|NFPT|silent|null", ignore_case=TRUE)) |
           !str_detect(characteristics, regex("functional|GH|PRL|ACTH|TSH|LH", ignore_case=TRUE)),
    group = ifelse(nfpt, "NFPT", "other_PitNET"),
    omics = "bulk_transcriptome_RNAseq"
  ) %>%
  select(cohort, geo_sample, platform, title, source, characteristics, tgroup, group, nfpt, omics)

master <- bind_rows(g1, g2, g3)
write_csv(master, file.path(OUT, "master_metadata.csv"))

cat("=== MASTER METADATA SUMMARY ===\n")
cat("Total samples:", nrow(master), "\n")
print(table(master$cohort, master$group, useNA="ifany"))
cat("\nNFPT subset counts:\n")
print(table(master$cohort, master$nfpt, useNA="ifany"))
cat("\nSaved ->", file.path(OUT, "master_metadata.csv"), "\n")
