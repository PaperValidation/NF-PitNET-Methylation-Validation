#!/usr/bin/env Rscript
suppressMessages({ library(readr); library(stringr) })
m <- read_csv("/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/samples/GSE207937_samples.csv", show_col_types=FALSE)
cat("=== unique subtype tokens ===\n")
print(table(str_extract(m$sample_characteristics_ch1, "(?<=subtypes:\\s)\\S+"), useNA="ifany"))
cat("\n=== raw characteristics head(15) ===\n")
for (i in 1:15) cat(i, ": ", m$sample_characteristics_ch1[i], "\n")
