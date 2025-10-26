#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GSVA); library(BiocParallel); library(GSEABase); library(readr); library(dplyr); library(tibble)
  library(msigdbr)
  library(AnnotationDbi); library(org.Mm.eg.db)   # for mapping if needed
})

#------------------------------
# GSVA (Hallmark gene sets)
#--------------------------------

proj   <- Sys.getenv("PROJECT", unset="~/glaucoma_rnaseq_pipeline")
indir  <- file.path(proj, "results/delivery")
outdir <- file.path(proj, "results/pathways")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Expression (genes x samples) — rows must be gene IDs
vst <- readr::read_tsv(file.path(indir,"counts_vst.tsv")) 
mat <- as.matrix(vst[,-1])
rownames(mat) <- vst$ENSEMBL

# strip version suffix if present (e.g., ENSMUSG00000000001.1 -> ENSMUSG00000000001)
rownames(mat) <- sub("\\.\\d+$", "", rownames(mat))

# Hallmark gene sets (mouse) via msigdbr
#   -> build ENSEMBL-based list for GSVA
m_df <- msigdbr(species = "mouse", collection = "H")

# Prefer ENSEMBL from msigdbr if present; else map SYMBOL -> ENSEMBL
id_col <- if ("ensembl_gene" %in% names(m_df)) "ensembl_gene" else "gene_symbol"

if (id_col == "gene_symbol") {
  sym2ens <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(m_df$gene_symbol),
                                   keytype = "SYMBOL", columns = c("ENSEMBL", "SYMBOL")) %>%
             distinct(SYMBOL, ENSEMBL) %>% filter(!is.na(ENSEMBL))
  m_df <- m_df %>% left_join(sym2ens, by = c("gene_symbol" = "SYMBOL")) %>%
         rename(ensembl_gene = ENSEMBL)
}

gsets_list <- m_df %>%
  transmute(gs_name, ensembl_gene = sub("\\.\\d+$","", ensembl_gene)) %>%  # just in case
  filter(!is.na(ensembl_gene)) %>%
  group_split(gs_name) %>%
  setNames(unique(m_df$gs_name)) %>%
  lapply(function(df) unique(df$ensembl_gene))

# write a GMT for reuse
gmt_path <- file.path(proj, "data/pathways/h.all.mouse.ensembl.gmt")
dir.create(dirname(gmt_path), showWarnings = FALSE, recursive = TRUE)
con <- file(gmt_path, "w"); on.exit(close(con), add = TRUE)
for (nm in names(gsets_list)) writeLines(paste(c(nm, "na", gsets_list[[nm]]), collapse = "\t"), con)

# Pick a parallel backend (or serial). On Linux/macOS, MulticoreParam is fine.
bp <- MulticoreParam(workers = 4)   
# On Windows use: bp <- SnowParam(workers = 4)

# GSVA: Use the new parameter object API
# 1. Create a parameter object using the gsvaParam() constructor
gsva_param <- GSVA::gsvaParam(exprData = mat,
                              geneSets = gsets_list,
                              kcdf = "Gaussian",
                              minSize = 10,
                              maxSize = 500,
                              maxDiff = TRUE)

# 2. Pass the parameter object to the gsva() function
gsva_es <- GSVA::gsva(gsva_param, BPPARAM = bp)

write_tsv(as.data.frame(gsva_es) %>% rownames_to_column("geneset"),
          file.path(outdir, "GSVA_hallmark_scores.tsv"))

#----------------------------------------------------------------
# quick differential GSVA (IOP vs control) by per-gene-set t-test
#----------------------------------------------------------------

# 1) Build a clean sample→condition table (TSV) from SRP metadata
srp_md <- readr::read_csv(file.path(proj, "data/metadata/SRP394552_metadata.csv"),
                          show_col_types = FALSE)

samples_tbl <- srp_md %>%
  mutate(
    sample = run_accession,
    condition = case_when(
      grepl("control", paste(treatment, experiment_title), ignore.case = TRUE) ~ "control",
      grepl("\\bIOP\\b", paste(treatment, experiment_title), ignore.case = TRUE) ~ "IOP",
      TRUE ~ NA_character_
    )
  ) %>%
  transmute(sample, condition) %>%
  distinct()

# (optional) write the two-column file for future runs
readr::write_tsv(samples_tbl, file.path(proj, "data/metadata/samples.txt"))

# 2) Use it (no error now)
meta <- samples_tbl

# 3) Ensure all GSVA columns are present, and order them
stopifnot(all(colnames(gsva_es) %in% meta$sample))
meta <- meta %>% filter(sample %in% colnames(gsva_es))
gsva_es <- gsva_es[, meta$sample, drop = FALSE]

# 4) Run per-geneset t-tests (IOP vs control)
grp <- setNames(meta$condition, meta$sample)
pvals <- apply(gsva_es, 1, function(x) {
  t.test(x[grp == "IOP"], x[grp == "control"])$p.value
})
padj <- p.adjust(pvals, method = "BH")

res <- tibble(geneset = rownames(gsva_es), pval = pvals, padj = padj) %>%
  arrange(padj)
readr::write_tsv(res, file.path(outdir, "GSVA_hallmark_ttests.tsv"))
