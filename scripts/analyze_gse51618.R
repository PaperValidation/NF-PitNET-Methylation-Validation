#!/usr/bin/env Rscript
# GSE51618: differential expression invasive vs non-invasive NFPA (microarray)
# Output: top DE genes, volcano-ready table, basic GO/KEGG via fgsea (if gene sets available)
suppressMessages({
  library(limma); library(data.table); library(dplyr); library(ggplot2)
})
BASE <- "/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
MAT  <- file.path(BASE, "GSE51618_series_matrix.txt.gz")
OUT  <- file.path(BASE, "results", "GSE51618")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## read matrix (skip ! header lines dynamically)
message("Reading matrix...")
# locate the table begin marker
con <- gzfile(MAT, "rt")
header_lines <- readLines(con, n = 200)
close(con)
begin_idx <- which(header_lines == "!series_matrix_table_begin")
data_start <- begin_idx + 1
message("Table begins at line ", data_start)
# column names are the first data line; ID_REF is col 1
mat <- fread(MAT, sep = "\t", header = TRUE, skip = data_start - 1, fill = TRUE, data.table = FALSE)
idcol <- colnames(mat)[1]
# ensure unique rownames
mat[[idcol]] <- make.unique(as.character(mat[[idcol]]))
rownames(mat) <- mat[[idcol]]
mat[[idcol]] <- NULL
mat <- as.matrix(mat)
mode(mat) <- "numeric"
message("Matrix dims: ", nrow(mat), " probes x ", ncol(mat), " samples")

## sample groups from master metadata
meta <- fread(file.path(BASE, "metadata", "master_metadata.csv"))
print(colnames(meta))
g1 <- meta %>% filter(cohort == "GSE51618")
g1$geo_sample <- gsub('"', '', g1$geo_sample)   # strip quoting from CSV
print(g1[, .(geo_sample, group)])
# align columns to matrix
samp <- intersect(colnames(mat), g1$geo_sample)
mat <- mat[, samp, drop = FALSE]
g1 <- g1 %>% filter(geo_sample %in% samp) %>% arrange(match(geo_sample, colnames(mat)))
group <- factor(as.character(g1$group))  # NFPA_invasive / NFPA_non-invasive / normal_pituitary
message("Groups: ", paste(table(group), collapse = " "))
grp_nm <- make.names(levels(group))
names(grp_nm) <- levels(group)
message("Group level names -> ", paste(grp_nm, collapse = ", "))

## design: invasive vs non-invasive (exclude normal for primary contrast)
tumor <- which(group %in% c("NFPA_invasive","NFPA_non-invasive"))
matT <- mat[, tumor, drop = FALSE]
grpT <- make.names(as.character(group[tumor]))
# clean: drop probes with NA / non-finite / zero variance
keep <- apply(matT, 1, function(x) all(is.finite(x)) && sd(x) > 0)
matT <- matT[keep, , drop = FALSE]
message("Probes after cleaning: ", nrow(matT))
design <- model.matrix(~ 0 + grpT)
colnames(design) <- make.names(unique(grpT))
contr <- makeContrasts(inv_vs_noninv = NFPA_invasive - NFPA_non.invasive, levels = design)
fit <- lmFit(matT, design)
fit <- contrasts.fit(fit, contr)
fit <- eBayes(fit, trend = TRUE)
tt <- topTable(fit, coef = "inv_vs_noninv", number = Inf, adjust = "fdr")
tt$ID_REF <- rownames(tt)
fwrite(tt, file.path(OUT, "DE_invasive_vs_noninvasive.csv"))

## annotation: map probes to symbols via GPL570 (if available locally) else GENE_SYMBOL column
# GSE51618 platform is GPL570 (HG-U133A). Try annotate from matrix !Platform data if present.
# Quick approach: read gene symbol from series matrix if a GENE_SYMBOL column exists.
message("Top DE genes (inv vs non-inv):")
print(head(tt[order(tt$adj.P.Val), c("ID_REF","logFC","AveExpr","P.Value","adj.P.Val")], 20))

## volcano
png(file.path(OUT, "volcano.png"), width = 900, height = 700)
with(tt, {
  cols <- ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "red", "grey")
  plot(logFC, -log10(adj.P.Val), col = cols, pch = 19, cex = 0.6,
       xlab = "log2 FC (invasive/non-invasive)", ylab = "-log10 FDR)",
       main = "GSE51618 NFPA: invasive vs non-invasive")
  abline(h = -log10(0.05), v = c(-1,1), lty = 2, col = "blue")
})
dev.off()
message("Saved DE table + volcano to ", OUT)
