#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport); library(readr); library(dplyr); library(tibble)
  library(DESeq2)
})

#-----------------------------------------------
#Export matrices (counts, normalized, VST, TPM)
#------------------------------------------------

proj <- Sys.getenv("PROJECT", unset = "~/glaucoma_rnaseq_pipeline")
salmon_dir <- file.path(proj, "results/salmon")
meta_file  <- file.path(proj, "data/metadata/SRP394552_metadata.csv")
tx2g_file  <- file.path(proj, "data/reference/tx2gene.csv")  # 2 cols: TXNAME, GENEID
outdir     <- file.path(proj, "results/delivery")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- metadata (wide) ----
meta <- readr::read_csv(meta_file, show_col_types = FALSE)
stopifnot(all(c("run_accession","treatment","experiment_title") %in% names(meta)))

# derive condition once, reproducibly
meta <- meta %>%
  mutate(
    sample = run_accession,
    condition = case_when(
      grepl("control", paste(treatment, experiment_title), ignore.case = TRUE) ~ "control",
      grepl("\\bIOP\\b", paste(treatment, experiment_title), ignore.case = TRUE) ~ "IOP",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition)) %>%
  mutate(condition = factor(condition, levels = c("control","IOP")))

# ---- files ----
stopifnot(all(file.exists(file.path(salmon_dir, meta$sample, "quant.sf"))))
files <- setNames(file.path(salmon_dir, meta$sample, "quant.sf"), meta$sample)

# ---- tx2gene (required for gene-level import) ----
tx2gene <- readr::read_csv(tx2g_file, show_col_types = FALSE)
stopifnot(ncol(tx2gene) >= 2)
names(tx2gene)[1:2] <- c("TXNAME","GENEID")
tx2gene <- tx2gene[,1:2]

# ---- tximport: gene-level (lengthScaledTPM is fine for DESeq2) ----
txi <- tximport(files, type = "salmon", txOut = FALSE, tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

# ---- DESeq2 object ----
coldata <- meta %>% select(sample, condition) %>% as.data.frame()
rownames(coldata) <- coldata$sample
stopifnot(identical(names(files), rownames(coldata)))  # strict alignment

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# keep nonzero rows, then size factors
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- estimateSizeFactors(dds)

# ---- matrices ----
cts_raw  <- round(counts(dds, normalized = FALSE))
cts_norm <- round(counts(dds, normalized = TRUE))
vst_mat  <- assay(vst(dds, blind = FALSE))
tpm_mat  <- txi$abundance  # gene-level TPM after tx2gene

# ---- save ----
write_tsv(as.data.frame(cts_raw)  %>% rownames_to_column("ENSEMBL"), file.path(outdir,"counts_raw.tsv"))
write_tsv(as.data.frame(cts_norm) %>% rownames_to_column("ENSEMBL"), file.path(outdir,"counts_norm.tsv"))
write_tsv(as.data.frame(vst_mat)  %>% rownames_to_column("ENSEMBL"), file.path(outdir,"counts_vst.tsv"))
write_tsv(as.data.frame(tpm_mat)  %>% rownames_to_column("ENSEMBL"), file.path(outdir,"tpm.tsv"))

saveRDS(dds, file.path(outdir, "dds.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_export.txt"))
