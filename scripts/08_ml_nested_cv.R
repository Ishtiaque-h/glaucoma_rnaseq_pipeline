#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(ranger)
  library(pROC)
  library(ggplot2)
})

# ---------------------------
# Fixed inputs
# ---------------------------
SAMPLES <- "data/metadata/samples.txt"

ACT_TABLES <- c(
  DoRothEA = "results/pathways/DoRothEA_TFactivity.tsv",
  GSVA     = "results/pathways/GSVA_hallmark_scores.tsv",
  PROGENy  = "results/pathways/PROGENy_scores.tsv"
)

OUTDIR <- "results/ml"

# ---------------------------
# Function to run one analysis
# ---------------------------
run_one <- function(NAME, SAMPLES, ACT_TABLE, OUTDIR) {
  cat("\n======================\nRunning:", NAME, "\n======================\n")

  dir.create(file.path(OUTDIR, NAME), recursive = TRUE, showWarnings = FALSE)
  OD <- file.path(OUTDIR, NAME)

  # ---------------------------
  # Load metadata
  # ---------------------------
  samples <- read.delim(SAMPLES, header = TRUE, check.names = FALSE) |>
    as_tibble()

  # If 'condition' exists, rename -> group
  if ("condition" %in% names(samples)) {
    samples <- samples |> dplyr::rename(group = condition)
  }

  stopifnot("sample" %in% names(samples))
  stopifnot("group"  %in% names(samples))

  # Force canonical labels IOP/control
  samples <- samples |>
    mutate(
      group_raw = group,
      group = case_when(
        str_detect(str_to_lower(group_raw), "iop|treated|treatment")     ~ "IOP",
        str_detect(str_to_lower(group_raw), "ctrl|control|untreated")    ~ "control",
        TRUE ~ NA_character_
      )
    )

  cat("Counts in samples after canonicalization:\n")
  print(table(samples$group, useNA = "ifany"))

  # ---------------------------
  # Load activity table, detect orientation
  # ---------------------------
  act_raw <- read.delim(ACT_TABLE, header = TRUE, check.names = FALSE)

  meta_ids       <- samples$sample
  first_col_vals <- act_raw[[1]]

  if ("sample" %in% names(act_raw) &&
      sum(act_raw$sample %in% meta_ids) > 0) {
    # Case 1: already has 'sample' column
    mat <- as_tibble(act_raw)

  } else if (sum(first_col_vals %in% meta_ids) > 0) {
    # Case 2: first column values are sample IDs (e.g. PROGENy)
    mat <- as_tibble(act_raw)
    colnames(mat)[1] <- "sample"

  } else {
    # Case 3: first column is features, rest are samples (e.g. DoRothEA / GSVA)
    feature_ids <- first_col_vals
    expr <- act_raw[, -1, drop = FALSE]
    rownames(expr) <- feature_ids

    mat <- t(as.matrix(expr)) |>
      as.data.frame() |>
      tibble::rownames_to_column("sample") |>
      as_tibble()
  }

  cat("Head of mat$sample (activity table samples):\n")
  print(head(mat$sample))

  # ---------------------------
  # Merge data
  # ---------------------------
  df <- samples |>
    select(sample, group) |>
    inner_join(mat, by = "sample")

  cat("Counts in df (after join with activity table):\n")
  print(table(df$group, useNA = "ifany"))

  cat("Samples in metadata but not in activity table:\n")
  print(setdiff(samples$sample, df$sample))

  cat("Samples in activity table but not in metadata:\n")
  print(setdiff(mat$sample, df$sample))

  # Drop rows with missing group
  df <- df |> filter(!is.na(group))

  # ---------------------------
  # Prepare X, y
  # ---------------------------
  # outcome: make 2-level factor with POSITIVE class first ("IOP")
  df$group <- factor(df$group, levels = c("IOP", "control"))
  stopifnot(all(levels(df$group) == c("IOP", "control")))

  # predictors
  X <- df |>
    select(-sample, -group)

  # drop all-NA columns if any
  X <- X[, colSums(!is.na(X)) > 0, drop = FALSE]

  # convert non-numeric to numeric where possible
  for (j in seq_len(ncol(X))) {
    if (!is.numeric(X[[j]])) {
      suppressWarnings({
        asnum <- as.numeric(X[[j]])
      })
      if (sum(is.na(asnum)) < length(asnum)) X[[j]] <- asnum
    }
  }

  # NZV removal
  nzv <- caret::nearZeroVar(X)
  if (length(nzv)) X <- X[, -nzv, drop = FALSE]

  # Drop columns that are mostly NA (optional but safer)
  na_frac <- colMeans(is.na(X))
  X <- X[, na_frac < 0.8, drop = FALSE]

  # complete cases
  keep <- stats::complete.cases(X)
  X <- X[keep, , drop = FALSE]
  y <- droplevels(df$group[keep])

  cat("Counts in y after complete.cases filtering:\n")
  print(table(y, useNA = "ifany"))

  taby <- table(y)
  if (length(taby) != 2 || any(taby < 2)) {
    stop("After filtering, not enough samples per class to train (need ≥2 in each).")
  }

  # ---------------------------
  # Resampling control (stratified, probabilistic, ROC-ready)
  # ---------------------------
  set.seed(123)
  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )

  # ---------------------------
  # GLMNET (logistic)
  # ---------------------------
  set.seed(42)
  fit_glm <- caret::train(
    x = X, y = y,
    method = "glmnet",
    family = "binomial",
    metric = "ROC",
    trControl = ctrl,
    preProcess = c("center", "scale"),
    tuneLength = 25
  )

  # ---------------------------
  # Ranger (RF)
  # ---------------------------
  set.seed(43)
  fit_rf <- caret::train(
    x = X, y = y,
    method = "ranger",
    metric = "ROC",
    trControl = ctrl,
    importance = "permutation",
    tuneGrid = expand.grid(
      mtry = max(1, floor(sqrt(ncol(X)))),
      splitrule = "gini",
      min.node.size = 1
    )
  )

  # ---------------------------
  # Metrics & ROC curves
  # ---------------------------
  summ_glm <- fit_glm$results %>% arrange(desc(ROC)) %>% slice(1)
  summ_rf  <- fit_rf$results  %>% arrange(desc(ROC))  %>% slice(1)

  write.csv(fit_glm$results,
            file.path(OD, paste0(NAME, "_glmnet_cv_results.csv")),
            row.names = FALSE)
  write.csv(fit_rf$results,
            file.path(OD, paste0(NAME, "_rf_cv_results.csv")),
            row.names = FALSE)

  cat("Best GLMNET ROC:", summ_glm$ROC, "  Best RF ROC:", summ_rf$ROC, "\n")

  # Probabilities from best models for ROC curve on pooled CV predictions
  pred_glm <- fit_glm$pred %>%
    dplyr::filter(lambda == fit_glm$bestTune$lambda,
                  alpha  == fit_glm$bestTune$alpha)

  pred_rf  <- fit_rf$pred

  roc_glm <- pROC::roc(response = pred_glm$obs,
                       predictor = pred_glm$IOP,
                       levels = c("control", "IOP"),
                       direction = "<")

  roc_rf  <- pROC::roc(response = pred_rf$obs,
                       predictor = pred_rf$IOP,
                       levels = c("control", "IOP"),
                       direction = "<")

  png(file.path(OD, paste0(NAME, "_glmnet_ROC.png")),
      width = 1200, height = 900, res = 160)
  plot(roc_glm,
       main = paste(NAME, "GLMNET ROC (AUC =", round(pROC::auc(roc_glm), 3), ")"))
  dev.off()

  png(file.path(OD, paste0(NAME, "_rf_ROC.png")),
      width = 1200, height = 900, res = 160)
  plot(roc_rf,
       main = paste(NAME, "RF ROC (AUC =", round(pROC::auc(roc_rf), 3), ")"))
  dev.off()

  # ---------------------------
  # Feature importance
  # ---------------------------
  vi_glm <- caret::varImp(fit_glm, scale = TRUE)$importance %>%
    rownames_to_column("feature") %>%
    arrange(desc(Overall))

  vi_rf <- caret::varImp(fit_rf, scale = TRUE)$importance %>%
    rownames_to_column("feature") %>%
    arrange(desc(Overall))

  write.csv(vi_glm,
            file.path(OD, paste0(NAME, "_glmnet_feature_importance.csv")),
            row.names = FALSE)
  write.csv(vi_rf,
            file.path(OD, paste0(NAME, "_rf_feature_importance.csv")),
            row.names = FALSE)

  plot_topk <- function(vi, title, out_png, k = 20) {
    p <- vi %>%
      slice_head(n = min(k, nrow(.))) %>%
      ggplot(aes(x = reorder(feature, Overall), y = Overall)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "Importance", title = title) +
      theme_minimal(base_size = 12)
    ggsave(out_png, p, width = 8, height = 8, dpi = 200)
  }

  plot_topk(vi_glm,
            paste(NAME, "GLMNET: Top Features"),
            file.path(OD, paste0(NAME, "_glmnet_feature_importance.png")))
  plot_topk(vi_rf,
            paste(NAME, "RF: Top Features"),
            file.path(OD, paste0(NAME, "_rf_feature_importance.png")))

  # Save fitted objects for reuse
  saveRDS(fit_glm, file.path(OD, paste0(NAME, "_glmnet_fit.rds")))
  saveRDS(fit_rf,  file.path(OD, paste0(NAME, "_rf_fit.rds")))

  message("Done: outputs in ", OD)
}

# ---------------------------
# Run for all three tables
# ---------------------------
for (nm in names(ACT_TABLES)) {
  run_one(nm, SAMPLES, ACT_TABLES[[nm]], OUTDIR)
}
