#!/usr/bin/env Rscript
BiocManager::install("progeny")
BiocManager::install("dorothea")
BiocManager::install("viper")
suppressPackageStartupMessages({
  library(progeny); library(dorothea); library(viper)
  library(readr); library(dplyr); library(tibble)
  library(org.Mm.eg.db); library(AnnotationDbi)
})

proj   <- Sys.getenv("PROJECT", unset="~/glaucoma_rnaseq_pipeline")
indir  <- file.path(proj, "results/delivery")
outdir <- file.path(proj, "results/pathways")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- Expression (VST is fine for activity inference) ---
vst  <- readr::read_tsv(file.path(indir, "counts_vst.tsv"))
expr <- as.matrix(vst[,-1])
rownames(expr) <- sub("\\.\\d+$","", vst$ENSEMBL)  # strip version if present

# --- Map ENSEMBL -> SYMBOL (Mouse), then collapse duplicates by mean ---
map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(rownames(expr)),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>% distinct(ENSEMBL, .keep_all = TRUE)

sym <- map$SYMBOL[match(rownames(expr), map$ENSEMBL)]
expr_sym <- expr[!is.na(sym) & sym != "", , drop = FALSE]
sym <- sym[!is.na(sym) & sym != ""]
rownames(expr_sym) <- sym

# collapse rows with the same SYMBOL (mean)
expr_sym <- rowsum(expr_sym, group = rownames(expr_sym)) / as.vector(table(rownames(expr_sym)))

# =========================
# PROGENy (mouse pathways)
# =========================
prog_scores <- progeny(
  expr_sym,
  scale    = TRUE,
  organism = "Mouse",
  top      = 500,
  perm     = 1000,
  z_scores = TRUE
)

write_tsv(
  as.data.frame(prog_scores) %>% rownames_to_column("pathway"),
  file.path(outdir, "PROGENy_scores.tsv")
)

# =========================
# DoRothEA (mouse regulons) + VIPER
# =========================
data(dorothea_mm, package = "dorothea")
regulons <- dorothea_mm %>%
  filter(confidence %in% c("A","B","C")) %>%
  dorothea::df2regulon()

tf_activity <- viper(
  expr_sym,
  regulons,
  verbose      = FALSE,
  eset.filter  = FALSE,
  scale        = TRUE
)

write_tsv(
  as.data.frame(tf_activity) %>% rownames_to_column("TF"),
  file.path(outdir, "DoRothEA_TFactivity.tsv")
)
