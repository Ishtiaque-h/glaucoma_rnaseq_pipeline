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

OUTDIR <- "results/ml/nested_cv"

# ---------------------------
# Function to run one analysis with nested CV
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

  if ("condition" %in% names(samples)) {
    samples <- samples |> dplyr::rename(group = condition)
  }

  stopifnot("sample" %in% names(samples))
  stopifnot("group"  %in% names(samples))

  samples <- samples |>
    mutate(
      group_raw = group,
      group = case_when(
        str_detect(str_to_lower(group_raw), "iop|treated|treatment")  ~ "IOP",
        str_detect(str_to_lower(group_raw), "ctrl|control|untreated") ~ "control",
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
    mat <- as_tibble(act_raw)
  } else if (sum(first_col_vals %in% meta_ids) > 0) {
    mat <- as_tibble(act_raw)
    colnames(mat)[1] <- "sample"
  } else {
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

  df <- df |> filter(!is.na(group))

  # ---------------------------
  # Prepare X, y (full data)
  # ---------------------------
  df$group <- factor(df$group, levels = c("IOP", "control"))
  stopifnot(all(levels(df$group) == c("IOP", "control")))

  X <- df |>
    select(-sample, -group)

  X <- X[, colSums(!is.na(X)) > 0, drop = FALSE]

  for (j in seq_len(ncol(X))) {
    if (!is.numeric(X[[j]])) {
      suppressWarnings({
        asnum <- as.numeric(X[[j]])
      })
      if (sum(is.na(asnum)) < length(asnum)) X[[j]] <- asnum
    }
  }

  nzv <- caret::nearZeroVar(X)
  if (length(nzv)) X <- X[, -nzv, drop = FALSE]

  na_frac <- colMeans(is.na(X))
  X <- X[, na_frac < 0.8, drop = FALSE]

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
  # Nested CV: outer folds
  # ---------------------------
  set.seed(2025)
  outer_folds <- caret::createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)

  outer_results <- tibble(
    fold  = integer(),
    model = character(),
    AUC   = numeric()
  )

  # store best hyperparameters per outer fold
  outer_best <- list()

  # clean empty tibble for predictions
  all_outer_pred <- tibble(
    fold     = integer(),
    model    = character(),
    obs      = factor(character(), levels = levels(y)),
    prob_IOP = numeric()
  )

  for (f in seq_along(outer_folds)) {
    cat("Outer fold:", f, "of", length(outer_folds), "\n")

    test_idx  <- outer_folds[[f]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)

    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx,  , drop = FALSE]
    y_test  <- y[test_idx]

    inner_ctrl <- caret::trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 5,
      classProbs = TRUE,
      summaryFunction = caret::twoClassSummary,
      savePredictions = "final"
    )

    # GLMNET (inner CV)
    set.seed(42 + f)
    fit_glm_inner <- caret::train(
      x = X_train, y = y_train,
      method = "glmnet",
      family = "binomial",
      metric = "ROC",
      trControl = inner_ctrl,
      preProcess = c("center", "scale"),
      tuneLength = 25
    )

    # RF (inner CV – currently single grid point)
    set.seed(43 + f)
    fit_rf_inner <- caret::train(
      x = X_train, y = y_train,
      method = "ranger",
      metric = "ROC",
      trControl = inner_ctrl,
      importance = "permutation",
      tuneGrid = expand.grid(
        mtry = max(1, floor(sqrt(ncol(X_train)))),
        splitrule = "gini",
        min.node.size = 1
      )
    )

    # ----- save bestTunes for transparency -----
    bt_glm <- fit_glm_inner$bestTune
    bt_glm$fold  <- f
    bt_glm$model <- "glmnet"
    bt_glm <- bt_glm[, c("fold", "model", setdiff(names(bt_glm), c("fold", "model")))]

    bt_rf <- fit_rf_inner$bestTune
    bt_rf$fold  <- f
    bt_rf$model <- "rf"
    bt_rf <- bt_rf[, c("fold", "model", setdiff(names(bt_rf), c("fold", "model")))]

    outer_best[[length(outer_best) + 1]] <- bt_glm
    outer_best[[length(outer_best) + 1]] <- bt_rf

    # Predictions on outer test fold
    prob_glm <- predict(fit_glm_inner, newdata = X_test, type = "prob")[, "IOP"]
    prob_rf  <- predict(fit_rf_inner,  newdata = X_test, type = "prob")[, "IOP"]

    roc_glm <- pROC::roc(response = y_test,
                         predictor = prob_glm,
                         levels = c("control", "IOP"),
                         direction = "<")
    roc_rf  <- pROC::roc(response = y_test,
                         predictor = prob_rf,
                         levels = c("control", "IOP"),
                         direction = "<")

    outer_results <- outer_results |>
      add_row(fold = f, model = "glmnet", AUC = as.numeric(pROC::auc(roc_glm))) |>
      add_row(fold = f, model = "rf",     AUC = as.numeric(pROC::auc(roc_rf)))

    all_outer_pred <- all_outer_pred |>
      bind_rows(
        tibble(
          fold      = f,
          model     = "glmnet",
          obs       = y_test,
          prob_IOP  = prob_glm
        ),
        tibble(
          fold      = f,
          model     = "rf",
          obs       = y_test,
          prob_IOP  = prob_rf
        )
      )
  }

  # Save nested CV outer AUCs
  write.csv(
    outer_results,
    file.path(OD, paste0(NAME, "_nested_outer_fold_auc.csv")),
    row.names = FALSE
  )

  outer_summary <- outer_results |>
    group_by(model) |>
    summarise(
      mean_AUC = mean(AUC),
      sd_AUC   = sd(AUC),
      .groups  = "drop"
    )

  write.csv(
    outer_summary,
    file.path(OD, paste0(NAME, "_nested_outer_mean_auc.csv")),
    row.names = FALSE
  )

  cat("Nested CV outer AUC summary for", NAME, ":\n")
  print(outer_summary)

  # ----- save outer best hyperparameters -----
  outer_best_tune <- bind_rows(outer_best)
  write.csv(
    outer_best_tune,
    file.path(OD, paste0(NAME, "_nested_outer_best_tune.csv")),
    row.names = FALSE
  )

  # ----- save all outer predictions -----
  write.csv(
    all_outer_pred,
    file.path(OD, paste0(NAME, "_nested_outer_predictions.csv")),
    row.names = FALSE
  )

  # ---------------------------
  # ROC plots from pooled outer predictions + CI
  # ---------------------------
  pred_glm_outer <- all_outer_pred |> filter(model == "glmnet")
  pred_rf_outer  <- all_outer_pred |> filter(model == "rf")

  roc_glm_outer <- pROC::roc(
    response  = pred_glm_outer$obs,
    predictor = pred_glm_outer$prob_IOP,
    levels    = c("control", "IOP"),
    direction = "<"
  )

  roc_rf_outer <- pROC::roc(
    response  = pred_rf_outer$obs,
    predictor = pred_rf_outer$prob_IOP,
    levels    = c("control", "IOP"),
    direction = "<"
  )

  png(file.path(OD, paste0(NAME, "_glmnet_nested_outer_ROC.png")),
      width = 1200, height = 900, res = 160)
  plot(roc_glm_outer,
       main = paste(NAME,
                    "GLMNET Nested Outer ROC (AUC =",
                    round(pROC::auc(roc_glm_outer), 3), ")"))
  dev.off()

  png(file.path(OD, paste0(NAME, "_rf_nested_outer_ROC.png")),
      width = 1200, height = 900, res = 160)
  plot(roc_rf_outer,
       main = paste(NAME,
                    "RF Nested Outer ROC (AUC =",
                    round(pROC::auc(roc_rf_outer), 3), ")"))
  dev.off()

  # ----- pooled outer AUC + 95% CI summary -----
  ci_glm_outer <- as.numeric(pROC::ci.auc(roc_glm_outer))
  ci_rf_outer  <- as.numeric(pROC::ci.auc(roc_rf_outer))

  pooled_auc_summary <- tibble(
    model       = c("glmnet", "rf"),
    auc_pooled  = c(as.numeric(pROC::auc(roc_glm_outer)),
                    as.numeric(pROC::auc(roc_rf_outer))),
    auc_ci_low  = c(ci_glm_outer[1], ci_rf_outer[1]),
    auc_ci_high = c(ci_glm_outer[3], ci_rf_outer[3])
  )

  write.csv(
    pooled_auc_summary,
    file.path(OD, paste0(NAME, "_nested_outer_pooled_auc_ci.csv")),
    row.names = FALSE
  )

  # ---------------------------
  # FINAL models on full data (for varImp, CV tables, etc.)
  # ---------------------------
  set.seed(123)
  ctrl_full <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )

  # GLMNET full-data
  set.seed(42)
  fit_glm_full <- caret::train(
    x = X, y = y,
    method = "glmnet",
    family = "binomial",
    metric = "ROC",
    trControl = ctrl_full,
    preProcess = c("center", "scale"),
    tuneLength = 25
  )

  # RF full-data
  set.seed(43)
  fit_rf_full <- caret::train(
    x = X, y = y,
    method = "ranger",
    metric = "ROC",
    trControl = ctrl_full,
    importance = "permutation",
    tuneGrid = expand.grid(
      mtry = max(1, floor(sqrt(ncol(X)))),
      splitrule = "gini",
      min.node.size = 1
    )
  )

  write.csv(
    fit_glm_full$results,
    file.path(OD, paste0(NAME, "_glmnet_cv_results_full_data.csv")),
    row.names = FALSE
  )
  write.csv(
    fit_rf_full$results,
    file.path(OD, paste0(NAME, "_rf_cv_results_full_data.csv")),
    row.names = FALSE
  )

  # ---------------------------
  # Feature importance from full-data models
  # ---------------------------
  vi_glm <- caret::varImp(fit_glm_full, scale = TRUE)$importance |>
    rownames_to_column("feature") |>
    arrange(desc(Overall))

  vi_rf <- caret::varImp(fit_rf_full, scale = TRUE)$importance |>
    rownames_to_column("feature") |>
    arrange(desc(Overall))

  write.csv(
    vi_glm,
    file.path(OD, paste0(NAME, "_glmnet_feature_importance.csv")),
    row.names = FALSE
  )
  write.csv(
    vi_rf,
    file.path(OD, paste0(NAME, "_rf_feature_importance.csv")),
    row.names = FALSE
  )

  plot_topk <- function(vi, title, out_png, k = 20) {
    top_vi <- vi %>%
      dplyr::slice_head(n = min(k, nrow(vi)))

    p <- top_vi %>%
      ggplot(aes(x = reorder(feature, Overall), y = Overall)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "Importance", title = title) +
      theme_minimal(base_size = 12)
    ggsave(out_png, p, width = 8, height = 8, dpi = 200)
  }

  plot_topk(
    vi_glm,
    paste(NAME, "GLMNET: Top Features (full-data model)"),
    file.path(OD, paste0(NAME, "_glmnet_feature_importance.png"))
  )
  plot_topk(
    vi_rf,
    paste(NAME, "RF: Top Features (full-data model)"),
    file.path(OD, paste0(NAME, "_rf_feature_importance.png"))
  )

  saveRDS(fit_glm_full, file.path(OD, paste0(NAME, "_glmnet_fit_full.rds")))
  saveRDS(fit_rf_full,  file.path(OD, paste0(NAME, "_rf_fit_full.rds")))

  message("Done: outputs in ", OD)
}

# ---------------------------
# Run for all three tables
# ---------------------------
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
for (nm in names(ACT_TABLES)) {
  run_one(nm, SAMPLES, ACT_TABLES[[nm]], OUTDIR)
}
