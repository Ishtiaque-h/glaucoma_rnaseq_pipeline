#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(effsize)       # Cliff's delta
})

# ---------------------------
# Config (edit paths if needed)
# ---------------------------
meta_file <- "data/metadata/samples.txt"
files <- list(
  DoRothEA = "results/pathways/DoRothEA_TFactivity.tsv",
  GSVA     = "results/pathways/GSVA_hallmark_scores.tsv",
  PROGENy  = "results/pathways/PROGENy_scores.tsv"
)
outdir <- "results/plots/activity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Helpers
# ---------------------------
read_meta <- function(f){
  m <- read.delim(f, header = TRUE, check.names = FALSE) %>% as_tibble()
  # standardize names if needed
  nms <- tolower(names(m))
  if(!all(c("sample","group") %in% nms)){
    # try to coerce first two columns
    if(ncol(m) >= 2){
      names(m)[1:2] <- c("sample","group")
    } else stop("samples.txt must have at least 2 columns: sample, group")
  } else {
    names(m) <- nms
  }
  m %>% mutate(group = as.factor(group))
}

# Read a matrix that may be oriented either way.
# Return a tibble with rows = samples, columns = features, and "sample" column.
read_activity <- function(path, samples){
  mat <- read.delim(path, header = TRUE, check.names = FALSE)
  # if first column looks like rownames, promote it
  if(!any(names(mat) == "feature") && !any(names(mat) == "gene")){
    if(all(!grepl("\\s", mat[[1]]))) {
      rn <- mat[[1]]
      mat <- mat[,-1, drop=FALSE]
      rownames(mat) <- rn
    }
  }
  # Decide orientation by matching names to sample IDs
  cn <- colnames(mat)
  rn <- rownames(mat)
  smp <- samples$sample

  if(!is.null(cn) && length(intersect(cn, smp)) > length(intersect(rownames(mat), smp))){
    # columns are samples -> transpose to rows=samples
    dat <- as_tibble(t(as.matrix(mat)), rownames = "sample")
    names(dat) <- make.names(names(dat), unique = TRUE)  # keep unique
  } else if(!is.null(rownames(mat)) && length(intersect(rownames(mat), smp)) > 0){
    dat <- as_tibble(mat, rownames = "sample")
  } else {
    stop(paste("Cannot infer sample orientation for", path))
  }

  # keep only samples of interest, order like metadata
  dat <- dat %>% filter(sample %in% smp)
  dat <- samples %>% select(sample) %>% left_join(dat, by="sample")
  dat
}

cohen_d <- function(x, y){
  m1 <- mean(x, na.rm = TRUE); m2 <- mean(y, na.rm = TRUE)
  s1 <- sd(x, na.rm = TRUE);   s2 <- sd(y, na.rm = TRUE)
  sp <- sqrt(((length(x)-1)*s1^2 + (length(y)-1)*s2^2) / (length(x)+length(y)-2))
  if(is.finite(sp) && sp > 0) (m1 - m2)/sp else NA_real_
}

# ---------------------------
# Main
# ---------------------------
meta <- read_meta(meta_file)

for(tag in names(files)){
  message("Processing: ", tag)
  X <- read_activity(files[[tag]], meta)

  long <- X %>%
    pivot_longer(-sample, names_to = "feature", values_to = "score") %>%
    left_join(meta, by = "sample")

  # stats per feature
  stats <- long %>%
    group_by(feature) %>%
    summarise(
      n_ctrl = sum(group == levels(group)[1]),
      n_iop  = sum(group == levels(group)[2]),
      mean_ctrl = mean(score[group == levels(group)[1]], na.rm = TRUE),
      mean_iop  = mean(score[group == levels(group)[2]], na.rm = TRUE),
      d = cohen_d(score[group == levels(group)[2]],
                  score[group == levels(group)[1]]),
      cliffs_delta = tryCatch(
        effsize::cliff.delta(score ~ group)$estimate, error = function(e) NA_real_),
      p_wilcox = tryCatch(
        wilcox.test(score ~ group)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(FDR = p.adjust(p_wilcox, method = "BH")) %>%
    arrange(FDR, desc(abs(d)))

  write.csv(stats, file.path(outdir, paste0(tag, "_effect_sizes.csv")), row.names = FALSE)

  # boxplots for top 20 by |d|
  top_feats <- stats %>% 
  slice_max(order_by = abs(d), n = 20, with_ties = FALSE) %>% 
  pull(feature)

  p <- long %>%
    filter(feature %in% top_feats) %>%
    mutate(feature = factor(feature, levels = top_feats)) %>%
    ggplot(aes(group, score, fill = group)) +
    geom_boxplot(outlier.size = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
    facet_wrap(~ feature, scales = "free_y") +
    scale_fill_manual(values = c("#989fb3","#4db6ac")) +
    labs(title = paste0(tag, " activity: top features by |Cohen's d|"),
         x = NULL, y = "Activity score") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))

  ggsave(file.path(outdir, paste0(tag, "_boxplots_top20.png")),
         p, width = 14, height = 10, dpi = 300)
}

message("All activity plots & stats written to: ", outdir)
