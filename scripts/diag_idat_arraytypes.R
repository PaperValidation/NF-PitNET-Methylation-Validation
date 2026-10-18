#!/usr/bin/env Rscript
source("C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/scripts/paths.R")
# Diagnostic: classify each IDAT array by probe count, map to SampleSheet rows.
suppressMessages(library(illuminaio))
suppressMessages(library(minfi))

BASE <- file.path(paths$DATASET_PATH_WSL, "Mendeley_Moller2025")
SS   <- file.path(BASE, "SampleSheet.xlsx")

# 1) read SampleSheet (first sheet)
# use openxlsx if available else readxl
if (requireNamespace("openxlsx", quietly=TRUE)) {
  sheet <- openxlsx::read.xlsx(SS, sheet=1)
} else if (requireNamespace("readxl", quietly=TRUE)) {
  sheet <- as.data.frame(readxl::read_excel(SS))
} else {
  stop("need openxlsx or readxl")
}
cat("SampleSheet rows:", nrow(sheet), "\n")
cat("SampleSheet cols:", paste(names(sheet), collapse=", "), "\n")

# 2) find idat files
idat_dir <- file.path(BASE, "Deidentified_IDATs")
grn <- list.files(idat_dir, pattern="_Grn.idat$", recursive=TRUE, full.names=TRUE)
cat("Total Grn idat found:", length(grn), "\n")

# 3) classify by probe count
classify <- function(n) {
  if (n < 700000)  return("450K")
  if (n < 1080000) return("EPICv1")
  return("EPICv2")
}
get_n <- function(b) {
  idat <- tryCatch(illuminaio::readIDAT(paste0(b, "_Grn.idat")), error=function(e) NULL)
  if (is.null(idat)) return(list(n=-1, chip="ERR"))
  n <- if (!is.null(idat$n)) idat$n else nrow(idat$Quants)
  chip <- if (!is.null(idat$ChipType)) idat$ChipType else NA_character_
  list(n=ifelse(is.null(n), -1, n), chip=chip)
}
types <- character(length(grn))
ns    <- integer(length(grn))
chips <- character(length(grn))
basenames <- character(length(grn))
for (i in seq_along(grn)) {
  b <- sub("_Grn.idat$", "", grn[i])
  basenames[i] <- b
  r <- get_n(b)
  ns[i] <- r$n
  chips[i] <- as.character(r$chip)
  types[i] <- if (r$n < 0) "ERR" else classify(r$n)
}
cat("\n=== Array type distribution (by probe count) ===\n")
print(table(types))
cat("probe-count range:", min(ns[ns>0]), "to", max(ns[ns>0]), "\n")

# 4) map basenames to SampleSheet Sentrix_ID + Sentrix_Position
if (all(c("Sentrix_ID","Sentrix_Position") %in% names(sheet))) {
  sheet$basename <- paste0(sheet$Sentrix_ID, "_", sheet$Sentrix_Position)
  # for each sheet row, find its array type (match by basename)
  sheet$type <- NA_character_
  for (j in seq_len(nrow(sheet))) {
    m <- which(basenames == sheet$basename[j])
    if (length(m)==1) sheet$type[j] <- types[m]
    else if (length(m)>1) sheet$type[j] <- paste(unique(types[m]), collapse="/")
  }
  cat("\n=== SampleSheet array types (per Sample_ID) ===\n")
  print(table(sheet$type, useNA="ifany"))
  # duplicate basenames?
  dup_base <- basenames[duplicated(basenames)]
  cat("\nDuplicate idat basenames in folder:", length(dup_base), "\n")
  # sheet basenames present multiple times
  dup_sheet <- sheet$basename[duplicated(sheet$basename)]
  cat("Duplicate basenames in SampleSheet:", length(dup_sheet), "\n")
  if (length(dup_sheet)>0) print(head(dup_sheet))
}
cat("\nDONE\n")
