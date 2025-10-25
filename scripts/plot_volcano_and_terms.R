library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(patchwork); library(stringr)

res <- read_csv("results/deseq2/deseq2_iop_vs_ctrl.csv", show_col_types = FALSE)
if (colnames(res)[1] == "") colnames(res)[1] <- "gene"
res <- res %>% mutate(
  neglog10p = -log10(pvalue + 1e-300),
  sig = padj < 0.05 & abs(log2FoldChange) > 1,
  direction = ifelse(log2FoldChange > 0, "Up (IOP)", "Down (IOP)")
)

# label top N by p-value in each direction
N <- 10
lab_up   <- res %>% filter(sig, log2FoldChange>0)  %>% arrange(padj) %>% slice_head(n=N)
lab_down <- res %>% filter(sig, log2FoldChange<0) %>% arrange(padj) %>% slice_head(n=N)
labs <- bind_rows(lab_up, lab_down)

p_volcano <- ggplot(res, aes(x=log2FoldChange, y=neglog10p, color=direction)) +
  geom_point(alpha=.5, size=1.2, na.rm=TRUE) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey60") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey60") +
  geom_text_repel(data=labs, aes(label=gene), size=3, max.overlaps=Inf, box.padding=.3, show.legend=FALSE) +
  scale_color_manual(values=c("Down (IOP)"="#2c7fb8","Up (IOP)"="#d95f0e")) +
  labs(title="IOP vs Control: Volcano plot", x="log2 fold-change (IOP/Control)", y="-log10(p-value)", color="Direction") +
  theme_bw()

# Read your enrichment (adjust paths if needed)
go_bp <- read.csv("results/pathways/ORA_GO_BP.csv")
kegg  <- read.csv("results/pathways/ORA_KEGG.csv")

fmt_terms <- function(tbl, title, n=10){
  top <- tbl %>% arrange(p.adjust) %>% slice_head(n=n) %>%
    mutate(term = paste0(gsub("_", " ", Description), "\nFDR=", signif(p.adjust,3)))
  ggplot(top, aes(x=reorder(term, p.adjust), y=-log10(p.adjust))) +
    geom_col() + coord_flip() + theme_bw() +
    labs(x=NULL, y="-log10(FDR)", title=title)
}

p_terms <- (fmt_terms(go_bp, "Top GO Biological Processes") | fmt_terms(kegg, "Top KEGG pathways"))

(p_volcano | p_terms) + plot_layout(widths=c(2,3))
ggsave("results/pathways/volcano_plus_terms.png", width=14, height=6, dpi=300)
