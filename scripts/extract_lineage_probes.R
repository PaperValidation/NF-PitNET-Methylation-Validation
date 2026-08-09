source("C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/scripts/paths.R")
suppressMessages({
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, quietly = TRUE)
  library(data.table, quietly = TRUE)
})

BASE <- paths$DATASET_PATH_WSL
OUT  <- file.path(BASE, "results", "GSE283928")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## EPIC annotation (probe-level)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- as.data.frame(anno)
anno_dt <- data.table(anno, keep.rownames = "Probe")
# gene name field
gcol <- grep("UCSC_RefGene_Name", colnames(anno_dt), value = TRUE)[1]
if (length(gcol) == 0) gcol <- "UCSC_RefGene_Name"
setnames(anno_dt, gcol, "Gene")

## Lineage marker gene sets for PitNET
## SF1-lineage (gonadotroph / silent gonadotroph / null cell): NR5A1(SF1), GATA3, GATA2, EMX2
## TPIT-lineage (corticotroph): TBX19(TPIT), POMC
## PIT1-lineage (somatotroph/lactotroph/thyrotroph): POU1F1(PIT1), GH1, PRL, TSHB
## Møller 2025 key genes (methylation signature anchors)
marker_sets <- list(
  SF1_lineage   = c("NR5A1", "GATA3", "GATA2", "EMX2", "FSHB", "LHB"),
  TPIT_lineage  = c("TBX19", "POMC", "PCSK1", "CRH"),
  PIT1_lineage  = c("POU1F1", "GH1", "PRL", "TSHB", "GHRHR"),
  Moller_key    = c("CCDC141", "ATP1A3", "ANK2", "CASC15", "HDC",
                    "NCKAP5", "GUCA1C", "ATP2C1", "GATA3", "NR5A1")
)

res <- list()
for (nm in names(marker_sets)) {
  genes <- marker_sets[[nm]]
  sub <- anno_dt[Gene %in% genes]
  # keep unique probes
  probes <- unique(sub$Probe)
  res[[nm]] <- data.table(set = nm, gene = sub$Gene, probe = sub$Probe)
  cat(sprintf("%-14s genes=%2d probes=%4d\n", nm, length(genes), length(probes)))
}

all_markers <- rbindlist(res)
fwrite(all_markers, file.path(OUT, "lineage_marker_probes.csv"))
cat("Saved lineage_marker_probes.csv with", nrow(all_markers), "probe-gene rows\n")
