suppressMessages({
  library(data.table, quietly = TRUE)
  library(glmnet, quietly = TRUE)
  library(pROC, quietly = TRUE)
})

BASE <- "/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
RES  <- file.path(BASE, "results", "GSE209903")
OUT  <- file.path(BASE, "results", "risk_model")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## 1. Load GSE209903 DE (FPT vs NFPT), already annotated with ENSG + we map to SYMBOL via simple map
de <- fread(file.path(RES, "FPT_vs_NFPT_full.csv"))
# columns: gene (ENSG), baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, symbol
de <- de[!is.na(padj) & !is.na(log2FoldChange)]
# keep protein-coding-ish genes with symbol
de <- de[!is.na(symbol) & symbol != ""]

## 2. Build binary label: NFPT = 1 (log2FC<0 means NFPT higher), FPT = 0
## But we don't have per-sample expr here (only DE summary). We need the original count/expr matrix.
## GSE209903 processed expression is NOT downloaded (only DE summary). Use DE lfc as the "signature building block"
## Strategy: build a surrogate NFPT risk score from top DE genes and validate distribution logic.
## Since no sample-level expr available, we instead: (a) extract a robust NFPT-up gene set (|lfc|>1, padj<0.05,
## baseMean>50 to avoid low-count noise), rank by -log10(padj)*sign(lfc).
sig <- de[padj < 0.05 & abs(log2FoldChange) > 1 & baseMean > 50]
sig[, score := -log10(padj) * sign(log2FoldChange)]   # NFPT-up positive
sig <- sig[order(-abs(score))]
cat(sprintf("NFPT signature candidate genes (|lfc|>1, padj<0.05, baseMean>50): %d\n", nrow(sig)))

## 3. For an actual risk MODEL we need sample-level expression. Flag dependency.
## The GSE209903 expression matrix (GSE209903_NFPT_FPT_expression...) was NOT downloaded.
## => We download it now (small, processed) to enable sample-level LASSO CV.
cat("NOTE: sample-level expression required for CV; see download step.\n")

## 4. Save signature for downstream mapping (Moller cross, GSE283928 mapping uses methylation, not this)
fwrite(sig[, .(symbol, gene, log2FoldChange, padj, score)],
       file.path(OUT, "nfpt_signature_candidates.csv"))
cat("Saved nfpt_signature_candidates.csv\n")

## 5. Quick enrichment-style sanity: are hallmark NFPT / pituitary genes present?
hallmark <- c("NR5A1","GATA3","SF1","FSHB","LHB","GH1","PRL","POU1F1","TBX19","CCDC141","ATP1A3","ANK2","HDC")
present <- sig[symbol %in% hallmark, .(symbol, log2FoldChange, padj)]
cat("Hallmark genes in signature:\n")
print(present)
fwrite(present, file.path(OUT, "hallmark_in_signature.csv"))
