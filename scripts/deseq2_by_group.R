#!/usr/bin/env Rscript
#BiocManager::install("pheatmap")
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

#project_dir  <- "glaucoma_rnaseq_pipeline"   # adjust if needed
meta_file    <- file.path("data/metadata", "metadata_deseq2.csv")
tx2gene_file <- file.path("data/reference", "tx2gene.csv")
salmon_dir   <- file.path("results/salmon")
outdir_root  <- file.path("results/deseq2")

dir.create(outdir_root, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Choose which comparison to run\
# EXAMPLES (you will need to uncomment each of these 3 examples, one at a time; so this script will need total 3 runs):
# -----------------------------

# For acute stress, uncomment this part
contrast_name    <- "Young_Acute_vs_Young_Ctrl"
groups_of_interest <- c("Young_Ctrl", "Young_Acute")      # baseline first, then treatment
platform_filter  <- "Illumina HiSeq 4000"                 # keep platform homogeneous

# For aging baseline:
#contrast_name    <- "Old_Ctrl_vs_Young_Ctrl"
#groups_of_interest <- c("Young_Ctrl", "Old_Ctrl")
#platform_filter  <- "Illumina HiSeq 4000"

# For repeated stress (NovaSeq only):
#contrast_name    <- "Young_Rep_IOP_vs_Young_Rep_Ctrl"
#groups_of_interest <- c("Young_Rep_Ctrl", "Young_Rep_IOP")
#platform_filter  <- "Illumina NovaSeq 6000"

outdir <- file.path(outdir_root, contrast_name)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Read metadata & subset
# -----------------------------
samples_all <- read.csv(meta_file, stringsAsFactors = FALSE)

samples <- samples_all %>%
  filter(Group_Label %in% groups_of_interest,
         instrument_model == platform_filter)

if (nrow(samples) == 0) {
  stop("No samples left after filtering. Check groups_of_interest and platform_filter.")
}

# Harmonize condition factor
samples$Group_Label <- factor(samples$Group_Label,
                              levels = groups_of_interest)

# Set rownames for colData
rownames(samples) <- samples$run_id

# -----------------------------
# 3. Link Salmon quant files
# -----------------------------
files <- file.path(salmon_dir, samples$run_id, "quant.sf")
names(files) <- samples$run_id
if (!all(file.exists(files))) {
  missing <- files[!file.exists(files)]
  stop("Missing Salmon quant files:\n", paste(missing, collapse = "\n"))
}

# -----------------------------
# 4. tximport: transcript -> gene
# -----------------------------
tx2gene <- read.csv(tx2gene_file, stringsAsFactors = FALSE)
colnames(tx2gene) <- tolower(colnames(tx2gene))
if (all(c("transcript_id","gene_id") %in% colnames(tx2gene))) {
  tx2gene <- tx2gene %>% rename(transcript = transcript_id, gene = gene_id)
} else if (!all(c("transcript","gene") %in% colnames(tx2gene))) {
  stop("tx2gene must have columns transcript/gene or transcript_id/gene_id")
}

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene[, c("transcript","gene")],
                countsFromAbundance = "lengthScaledTPM",
		ignoreTxVersion = TRUE)

# -----------------------------
# 5. Build DESeqDataSet (per run)
# -----------------------------
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ Group_Label)

# -----------------------------
# 6. Collapse technical replicates by biosample
# -----------------------------
dds_collapsed <- collapseReplicates(
  dds,
  groupby = dds$biosample_id,
  run     = dds$run_id
)

# Build new colData: one row per biosample
bio_ids <- colnames(dds_collapsed)

meta_collapsed <- samples %>%
  group_by(biosample_id) %>%
  summarise(
    Group_Label     = unique(Group_Label),
    instrument_model = unique(instrument_model)
  ) %>%
  ungroup()

stopifnot(setequal(bio_ids, meta_collapsed$biosample_id))
meta_collapsed <- meta_collapsed[match(bio_ids, meta_collapsed$biosample_id), ]
rownames(meta_collapsed) <- meta_collapsed$biosample_id

colData(dds_collapsed) <- S4Vectors::DataFrame(meta_collapsed)

# Ensure Group_Label factor levels are correct (baseline first)
dds_collapsed$Group_Label <- factor(dds_collapsed$Group_Label,
                                    levels = groups_of_interest)

# -----------------------------
# 7. Filter lowly expressed genes
# -----------------------------
keep <- rowSums(counts(dds_collapsed) >= 10) >= 2
dds_collapsed <- dds_collapsed[keep, ]

# -----------------------------
# 8. Run DESeq2
# -----------------------------
dds_collapsed <- DESeq(dds_collapsed)

# We want: groups_of_interest[2] vs groups_of_interest[1]
res <- results(dds_collapsed,
               contrast = c("Group_Label",
                            groups_of_interest[2],
                            groups_of_interest[1]))

coef_name <- paste0("Group_Label_", groups_of_interest[2],
                    "_vs_", groups_of_interest[1])

res_shrunk <- tryCatch(
  lfcShrink(dds_collapsed,
            coef = coef_name,
            type = "apeglm"),
  error = function(e) {
    message("apeglm not available; using normal shrinkage")
    lfcShrink(dds_collapsed,
              coef = coef_name,
              type = "normal")
  }
)

# -----------------------------
# 9. Annotate with gene symbols
# -----------------------------
res_df <- as.data.frame(res_shrunk)
res_df$ensembl <- rownames(res_df)

res_df$symbol <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = res_df$ensembl,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

res_df <- res_df[order(res_df$padj), ]

write.csv(res_df,
          file = file.path(outdir, paste0("deseq2_", contrast_name, "_annotated.csv")),
          row.names = FALSE)

sig <- subset(res_df, !is.na(padj) & padj < 0.05)
write.csv(sig,
          file = file.path(outdir, paste0("deseq2_", contrast_name, "_sig_padj0.05.csv")),
          row.names = FALSE)

# -----------------------------
# 10. PCA on collapsed samples
# -----------------------------
vsd <- vst(dds_collapsed, blind = FALSE)

p_pca <- plotPCA(vsd, intgroup = "Group_Label") +
  ggtitle(paste0("PCA (vst) – ", contrast_name)) +
  theme_bw()

ggsave(filename = file.path(outdir, paste0("pca_", contrast_name, ".png")),
       plot = p_pca, width = 6, height = 5, dpi = 150)

# -----------------------------
# 11. MA plot (shrunken LFC)
# -----------------------------
png(file.path(outdir, paste0("MA_plot_", contrast_name, ".png")),
    width = 1200, height = 900, res = 150)
plotMA(res_shrunk, ylim = c(-5,5),
       main = paste0("MA plot – ", contrast_name))
dev.off()

# -----------------------------
# 12. Volcano plot
# -----------------------------
volcano_df <- res_df %>%
  mutate(
    padj_safe = ifelse(is.na(padj) | padj == 0, NA, padj),
    negLog10Padj = -log10(padj_safe),
    sig = ifelse(!is.na(padj_safe) &
                   padj_safe < 0.05 &
                   abs(log2FoldChange) > 1,
                 "padj<0.05 & |LFC|>1", "ns")
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("ns" = "grey70",
                                "padj<0.05 & |LFC|>1" = "red")) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = paste0("Volcano: ", contrast_name),
       x = "log2 fold change (shrunken)",
       y = "-log10 adjusted p-value",
       color = "Significance") +
  theme_bw()

ggsave(filename = file.path(outdir, paste0("volcano_", contrast_name, ".png")),
       plot = p_volcano, width = 6, height = 5, dpi = 150)

# -----------------------------
# 13. Heatmap of top DE genes
# -----------------------------
top_n <- 30
top_genes <- head(rownames(res_df[order(res_df$padj), ]), top_n)

mat <- assay(vsd)[top_genes, ]

annotation_col <- data.frame(
  Group_Label = colData(vsd)$Group_Label
)
rownames(annotation_col) <- colnames(vsd)

pheatmap(mat,
         scale = "row",  # z-score each gene
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = paste0("Top DE genes – ", contrast_name),
         filename = file.path(outdir, paste0("heatmap_top_genes_", contrast_name, ".png")),
         width = 7, height = 7)

message("DESeq2 analysis completed for: ", contrast_name,
        "  (results in: ", outdir, ")")
