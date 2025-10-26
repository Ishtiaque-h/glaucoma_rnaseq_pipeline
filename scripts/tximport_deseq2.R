# install.packages("BiocManager")
#BiocManager::install(c("tximport","DESeq2","readr"))
#BiocManager::install(c("apeglm"))
library(apeglm)
library(tximport)
library(readr)
library(DESeq2)

# 1) metadata
samples <- read.csv("data/metadata/SRP394552_metadata.csv", stringsAsFactors=FALSE)
rownames(samples) <- samples$run_accession

# 2) point to Salmon quant files
files <- file.path("results/salmon", samples$run_accession, "quant.sf")
stopifnot(all(file.exists(files)))
names(files) <- samples$run_accession

# 3) tx2gene map
tx2gene <- read.csv("data/reference/tx2gene.csv", stringsAsFactors=FALSE)

# 4) import
txi <- tximport(files, type="salmon", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM",ignoreTxVersion = TRUE)

# 5) choose a clean first contrast: IOP vs CTRL (same age/platform if you want to be strict)
samples$condition <- ifelse(
  grepl("control", paste(samples$treatment, samples$experiment_title), ignore.case = TRUE), "control",
  ifelse(grepl("\\bIOP\\b", paste(samples$treatment, samples$experiment_title), ignore.case = TRUE), "IOP", NA)
)
samples$condition <- factor(samples$condition, levels = c("control","IOP"))
write.csv(samples, "SRP394552_metadata.csv", row.names = FALSE)
dds <- DESeqDataSetFromTximport(txi, colData=samples, design = ~ condition)

# (optional) filter low counts
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

# 6) run DE
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","IOP","control"))
res <- lfcShrink(dds, coef="condition_IOP_vs_control", type="apeglm")

# 7) write outputs
dir.create("results/deseq2", showWarnings=FALSE, recursive=TRUE)
write.csv(as.data.frame(res), "results/deseq2/deseq2_iop_vs_ctrl.csv")

# Quick QC: PCA
vsd <- vst(dds, blind=FALSE)
p <- plotPCA(vsd, intgroup="condition")
ggplot2::ggsave("results/deseq2/pca_condition.png", p, width=6, height=5, dpi=150)
