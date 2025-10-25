# install.packages("tidyverse")
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler","org.Mm.eg.db","enrichplot", "ReactomePA","msigdbr","DOSE","tidyverse"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)   # Mouse
  library(AnnotationDbi)
  library(enrichplot)
  library(ReactomePA)
  library(msigdbr)
  library(DOSE)
})

# --------------------------
# Inputs
# --------------------------
in_csv <- "results/deseq2/deseq2_iop_vs_ctrl.csv"   # adjust if different
outdir <- "results/pathways"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --------------------------
# Read DE results
# Expected columns: log2FoldChange, pvalue, padj
# Gene IDs in 1st column with no name and id type is "SYMBOL"
# --------------------------
res <- read.csv(in_csv)
colnames(res)[1] <- "gene"
gene_raw <- colnames(res)[1]
id_type <- "ENSEMBL"

# --------------------------
# Map to ENTREZ for enrichment
# --------------------------
map <- AnnotationDbi::select(org.Mm.eg.db,
                             keys = unique(res[[gene_raw]]),
                             keytype = id_type,
                             columns = c("ENTREZID","SYMBOL","ENSEMBL"))
map <- map %>% distinct(!!sym(id_type), .keep_all = TRUE)

res2 <- res %>%
  left_join(map, by = setNames(id_type, gene_raw))

# Keep a clean key
res2 <- res2 %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>%
  filter(!is.na(ENTREZID))

# --------------------------
# Define gene sets
# --------------------------
sig <- res2 %>% filter(padj < 0.05, abs(log2FoldChange) > 1)

# Ranked vector for GSEA (ENTREZID -> ranking metric)
# continuous score + tiny tie-breaker
rank_metric <- with(res2,
  sign(log2FoldChange) * (-log10(pvalue + 1e-300)) +
  1e-6 * rank(abs(log2FoldChange), ties.method = "random")
)

# one value per ENTREZID
rnk_df <- dplyr::tibble(ENTREZID = res2$ENTREZID, score = rank_metric) %>%
  dplyr::filter(!is.na(ENTREZID), is.finite(score)) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

rnk <- rnk_df$score
names(rnk) <- rnk_df$ENTREZID
rnk <- sort(rnk, decreasing = TRUE)

# --------------------------
# ORA: GO (BP, MF, CC)
# --------------------------
ego_bp <- enrichGO(gene = sig$ENTREZID, OrgDb = org.Mm.eg.db, keyType="ENTREZID",
                   ont="BP", pAdjustMethod="BH", readable=TRUE)
ego_mf <- enrichGO(gene = sig$ENTREZID, OrgDb = org.Mm.eg.db, keyType="ENTREZID",
                   ont="MF", pAdjustMethod="BH", readable=TRUE)
ego_cc <- enrichGO(gene = sig$ENTREZID, OrgDb = org.Mm.eg.db, keyType="ENTREZID",
                   ont="CC", pAdjustMethod="BH", readable=TRUE)

# --------------------------
# ORA: KEGG (mouse = 'mmu')
# --------------------------
ekegg <- enrichKEGG(gene = sig$ENTREZID, organism = "mmu", pAdjustMethod = "BH")
ekegg <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

# --------------------------
# ORA: Reactome
# --------------------------
ereact <- enrichPathway(gene = sig$ENTREZID, organism = "mouse", pAdjustMethod="BH", readable=TRUE)

# --------------------------
# GSEA: GO BP & KEGG
# --------------------------
gse_bp <- gseGO(geneList = rnk, OrgDb = org.Mm.eg.db, ont="BP",
                keyType="ENTREZID", pAdjustMethod="BH", verbose=FALSE)
gse_kegg <- gseKEGG(geneList = rnk, organism="mmu", verbose=FALSE)

# --------------------------
# Save tables
# --------------------------
write.csv(as.data.frame(ego_bp),  file.path(outdir, "ORA_GO_BP.csv"), row.names=FALSE)
write.csv(as.data.frame(ego_mf),  file.path(outdir, "ORA_GO_MF.csv"), row.names=FALSE)
write.csv(as.data.frame(ego_cc),  file.path(outdir, "ORA_GO_CC.csv"), row.names=FALSE)
write.csv(as.data.frame(ekegg),   file.path(outdir, "ORA_KEGG.csv"),  row.names=FALSE)
write.csv(as.data.frame(ereact),  file.path(outdir, "ORA_Reactome.csv"), row.names=FALSE)

write.csv(as.data.frame(gse_bp),   file.path(outdir, "GSEA_GO_BP.csv"),   row.names=FALSE)
write.csv(as.data.frame(gse_kegg), file.path(outdir, "GSEA_KEGG.csv"),    row.names=FALSE)

# --------------------------
# Save plots (top terms)
# --------------------------
plot_save <- function(p, name, w=8, h=6) {
  ggsave(filename = file.path(outdir, name), plot = p, width = w, height = h, dpi = 150)
}

if (nrow(as.data.frame(ego_bp)) > 0) plot_save(dotplot(ego_bp, showCategory=20), "ORA_GO_BP_dotplot.png")
if (nrow(as.data.frame(ekegg))  > 0) plot_save(dotplot(ekegg,  showCategory=20), "ORA_KEGG_dotplot.png")
if (nrow(as.data.frame(ereact)) > 0) plot_save(dotplot(ereact, showCategory=20), "ORA_Reactome_dotplot.png")

if (nrow(as.data.frame(gse_bp))   > 0) plot_save(gseaplot2(gse_bp,   geneSetID = gse_bp@result$ID[1]), "GSEA_GO_BP_leading_edge.png", w=9)
if (nrow(as.data.frame(gse_kegg)) > 0) plot_save(gseaplot2(gse_kegg, geneSetID = gse_kegg@result$ID[1]), "GSEA_KEGG_leading_edge.png", w=9)

message("Done. Results in: ", outdir)
