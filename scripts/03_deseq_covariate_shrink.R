suppressPackageStartupMessages({
  library(DESeq2); library(apeglm); library(dplyr); library(tibble)
  library(AnnotationDbi); library(org.Mm.eg.db); library(readr)
})

#----------------------------------------------
# DE with       covariate + LFC shrinkage
#----------------------------------------------
proj   <- Sys.getenv("PROJECT", unset = "~/glaucoma_rnaseq_pipeline")
outdir <- file.path(proj, "results/delivery")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dds <- readRDS(file.path(outdir, "dds.rds"))

#  add a covariate
colData(dds)$age_group <- factor(readr::read_csv(file.path(proj,"data/metadata/SRP394552_metadata.csv"))$age)
design(dds) <- ~ age_group + condition

#design(dds) <- ~ condition
dds <- DESeq(dds)

coef_name <- "condition_IOP_vs_control"
res    <- results(dds, contrast = c("condition","IOP","control"))
res_sh <- lfcShrink(dds, coef = coef_name, type = "apeglm")

res_df <- as.data.frame(res_sh) %>% rownames_to_column("ENSEMBL")

anno <- AnnotationDbi::select(org.Mm.eg.db,
  keys = res_df$ENSEMBL, keytype = "ENSEMBL",
  columns = c("SYMBOL","GENENAME","ENTREZID"))
anno <- dplyr::distinct(anno, ENSEMBL, .keep_all = TRUE)

res_anno <- res_df %>%
  left_join(anno, by="ENSEMBL") %>%
  arrange(padj, dplyr::desc(abs(log2FoldChange)))

write_csv(res_anno, file.path(outdir,"DE_full.csv"))
res_sig <- res_anno %>% filter(!is.na(padj), padj < 0.05)
write_csv(res_sig, file.path(outdir,"DE_sig.csv"))
write_csv(res_sig %>% filter(abs(log2FoldChange) >= 1), file.path(outdir,"DE_sig_absLFC1.csv"))

saveRDS(res_anno, file.path(outdir,"DE_full.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_deseq.txt"))
