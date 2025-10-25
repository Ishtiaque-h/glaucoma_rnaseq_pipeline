ego_bp <- read.csv("results/pathways/ORA_GO_BP.csv")#install.packages(c("igraph","ggplot2"))
# Bioc:
#if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#BiocManager::install(c("enrichplot","clusterProfiler","org.Mm.eg.db"))

library(clusterProfiler); library(enrichplot); library(ggplot2)

ego_bp <- readRDS("results/pathways/ORA_GO_BP.rds")

# (optional but useful) reduce redundancy before similarities
ego_bp <- clusterProfiler::simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

# compute term–term similarity and plot
ego_bp <- enrichplot::pairwise_termsim(ego_bp)

# if too few terms for a map, skip gracefully
if (nrow(as.data.frame(ego_bp)) >= 3) {
  p1 <- enrichplot::emapplot(ego_bp, showCategory = 30) +
        theme(text = element_text(size = 10))
  ggsave("results/pathways/emap_GO_BP.png", p1, width = 10, height = 8, dpi = 300)
} else {
  message("Not enough enriched terms for emapplot; try lowering p.adjust cutoff or increasing showCategory.")
}

