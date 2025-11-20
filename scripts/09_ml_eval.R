#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(pROC)
  library(ggplot2)
})

SAMPLES <- "data/metadata/samples.txt"
ACT_TABLES <- c(
  DoRothEA = "results/pathways/DoRothEA_TFactivity.tsv",
  GSVA     = "results/pathways/GSVA_hallmark_scores.tsv",
  PROGENy  = "results/pathways/PROGENy_scores.tsv"
)
INDIR  <- "results/ml"              # where your *_fit.rds live
OUTDIR <- "results/ml_eval"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers (identical merge/orientation logic to your script) ----------
load_samples <- function(path) {
  samples <- read.delim(path, header = TRUE, check.names = FALSE) |> as_tibble()
  if ("condition" %in% names(samples)) samples <- samples |> dplyr::rename(group = condition)
  stopifnot(all(c("sample","group") %in% names(samples)))
  samples |>
    mutate(group_raw = group,
           group = case_when(
             str_detect(str_to_lower(group_raw), "iop|treated|treatment") ~ "IOP",
             str_detect(str_to_lower(group_raw), "ctrl|control|untreated") ~ "control",
             TRUE ~ NA_character_
           )) |>
    filter(!is.na(group))
}

load_activity <- function(tab_path, samples) {
  act_raw <- read.delim(tab_path, header = TRUE, check.names = FALSE)
  meta_ids <- samples$sample
  first_col_vals <- act_raw[[1]]

  if ("sample" %in% names(act_raw) && sum(act_raw$sample %in% meta_ids) > 0) {
    mat <- as_tibble(act_raw)
  } else if (sum(first_col_vals %in% meta_ids) > 0) {
    mat <- as_tibble(act_raw); colnames(mat)[1] <- "sample"
  } else {
    feature_ids <- first_col_vals
    expr <- act_raw[, -1, drop = FALSE]
    rownames(expr) <- feature_ids
    mat <- t(as.matrix(expr)) |>
      as.data.frame() |>
      tibble::rownames_to_column("sample") |>
      as_tibble()
  }

  df <- samples |> select(sample, group) |> inner_join(mat, by = "sample")
  # predictors like your training:
  X <- df |> select(-sample, -group)
  X <- X[, colSums(!is.na(X)) > 0, drop = FALSE]
  for (j in seq_len(ncol(X))) {
    if (!is.numeric(X[[j]])) {
      suppressWarnings(asnum <- as.numeric(X[[j]]))
      if (sum(is.na(asnum)) < length(asnum)) X[[j]] <- asnum
    }
  }
  nzv <- caret::nearZeroVar(X); if (length(nzv)) X <- X[, -nzv, drop = FALSE]
  na_frac <- colMeans(is.na(X)); X <- X[, na_frac < 0.8, drop = FALSE]
  keep <- stats::complete.cases(X)
  list(df = df[keep, c("sample","group")], X = X[keep, , drop = FALSE])
}

best_cv_preds <- function(fit, positive = "IOP") {
  preds <- fit$pred

  # 1) If tuning params exist in preds, keep only bestTune rows
  if (!is.null(fit$bestTune)) {
    tune_cols <- intersect(names(fit$bestTune), names(preds))
    if (length(tune_cols)) {
      # join-style filter using all tune columns
      key <- fit$bestTune[1, tune_cols, drop = FALSE]
      for (cc in tune_cols) {
        if (is.numeric(preds[[cc]]) && is.numeric(key[[cc]])) {
          preds <- dplyr::filter(preds, abs(.data[[cc]] - key[[cc]]) < 1e-9)
        } else {
          preds <- dplyr::filter(preds, .data[[cc]] == key[[cc]])
        }
      }
    }
  }
  
  preds <- preds |>
  dplyr::arrange(Resample, rowIndex) |>
  dplyr::distinct(Resample, rowIndex, .keep_all = TRUE)


  # 2) Find the probability column for the POSITIVE class robustly
  prob_col <- if (positive %in% names(preds)) {
    positive
  } else {
    # try case-insensitive / punctuation variants
    cand <- grep(paste0("^\\W*", positive, "\\W*$"), names(preds), ignore.case = TRUE, value = TRUE)
    if (length(cand)) cand[1] else {
      stop("Could not find probability column for class '", positive, "'. Columns are: ",
           paste(names(preds), collapse = ", "))
    }
  }

  preds |>
    transmute(rowIndex, Resample,
              obs = .data[["obs"]],
              prob = .data[[prob_col]])
}

eval_from_preds <- function(preds, positive = "IOP", name = "MODEL", outdir = ".") {
  # Clean & validate
  preds <- preds |>
    dplyr::filter(!is.na(obs), !is.na(prob)) |>
    dplyr::mutate(
      obs = factor(as.character(obs), levels = c("control", positive))
    ) |>
    dplyr::filter(!is.na(obs))

  cls <- levels(droplevels(preds$obs))
  if (length(unique(preds$obs)) < 2) {
    stop("After cleaning, predictions contain a single class (",
         paste(cls, collapse = ","), ").")
  }
  if (!all(c("control", positive) %in% cls)) {
    stop("Observed classes must include both 'control' and '", positive,
         "'. Got: ", paste(cls, collapse = ","))
  }

  # ROC
  roc_obj <- pROC::roc(
    response  = preds$obs,
    predictor = preds$prob,
    levels    = c("control", positive),
    direction = "<",
    quiet     = TRUE
  )

  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_auc  <- as.numeric(pROC::ci.auc(roc_obj))

  # Threshold + derived metrics (Youden's J)
  thr <- pROC::coords(
    roc_obj, "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "ppv", "npv")
  )

  # Make sure thr is a data frame so we can index safely
  thr_df <- as.data.frame(thr)

  th_val <- as.numeric(thr_df$threshold[1])

  pred_label <- factor(
    ifelse(preds$prob >= th_val, positive, "control"),
    levels = c("control", positive)
  )
  stopifnot(length(pred_label) == length(preds$obs))

  cm <- caret::confusionMatrix(pred_label, preds$obs, positive = positive)

  # Summary tibble
  tibble::tibble(
    model       = name,
    auc         = round(auc_val, 3),
    auc_ci_low  = round(ci_auc[1], 3),
    auc_ci_high = round(ci_auc[3], 3),
    threshold   = as.numeric(thr_df$threshold[1]),
    sens        = as.numeric(thr_df$sensitivity[1]),
    spec        = as.numeric(thr_df$specificity[1]),
    ppv         = as.numeric(thr_df$ppv[1]),
    npv         = as.numeric(thr_df$npv[1]),
    accuracy    = cm$overall[["Accuracy"]],
    kappa       = cm$overall[["Kappa"]]
  ) |>
    write.csv(file.path(outdir, paste0(name, "_summary.csv")), row.names = FALSE)

  # Save ROC figure
  png(file.path(outdir, paste0(name, "_ROC_eval.png")),
      width = 1200, height = 900, res = 160)
  plot(
    roc_obj,
    main = paste0(
      name, " ROC (AUC=", round(auc_val, 3),
      ", 95% CI ", round(ci_auc[1], 3), "–", round(ci_auc[3], 3), ")"
    )
  )
  abline(a = 0, b = 1, lty = 2, col = "grey60")
  dev.off()

  # Save confusion matrix
  capture.output(cm) |>
    writeLines(con = file.path(outdir, paste0(name, "_confusion.txt")))

  invisible(list(roc = roc_obj, cm = cm, thr = thr_df))
}

rf_stability <- function(X, y, bestTune, runs=50, seed0=100, outdir=".", label="RF_stability") {
  set.seed(seed0)
  res <- vector("list", runs)
  ctrl_none <- caret::trainControl(method="none", classProbs=TRUE)
  tg <- bestTune
  for (i in seq_len(runs)) {
    set.seed(seed0 + i)
    fit <- caret::train(x = X, y = y,
                        method = "ranger",
                        trControl = ctrl_none,
                        importance = "permutation",
                        tuneGrid = tg,
                        metric = "Accuracy")
    vi <- caret::varImp(fit, scale=TRUE)$importance |>
      rownames_to_column("feature") |>
      arrange(desc(Overall)) |>
      mutate(run=i, rank=row_number())
    res[[i]] <- vi
  }
  all_vi <- bind_rows(res)
  freq <- all_vi |> group_by(feature) |>
    summarise(mean_rank = mean(rank),
              appear_top10 = mean(rank <= 10),
              .groups="drop") |>
    arrange(mean_rank)
  write.csv(freq, file.path(outdir, paste0(label, "_top10_frequency.csv")), row.names = FALSE)

  p <- freq |> slice_head(n=20) |>
    ggplot(aes(x=reorder(feature, -appear_top10), y=appear_top10)) +
    geom_col() + coord_flip() +
    labs(x=NULL, y="Frequency in top-10 over runs",
         title=paste0(label, " – feature stability (n=", runs, ")")) +
    theme_minimal(base_size=12)
  ggsave(file.path(outdir, paste0(label, "_top10_frequency.png")), p, width=8, height=6, dpi=200)
}

# ---------- main per-table evaluation ----------
eval_table <- function(NAME) {
  message("\n==== Evaluating ", NAME, " ====")
  table_path <- ACT_TABLES[[NAME]]
  model_dir  <- file.path(INDIR, NAME)
  outdir     <- file.path(OUTDIR, NAME); dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

  # Data (to map rowIndex -> sample and for stability refits)
  samples <- load_samples(SAMPLES)
  dat <- load_activity(table_path, samples)
  df_keep <- dat$df
  X <- dat$X
  y <- factor(df_keep$group, levels = c("IOP","control"))

  # Load fits
  fit_glm <- readRDS(file.path(model_dir, paste0(NAME, "_glmnet_fit.rds")))
  fit_rf  <- readRDS(file.path(model_dir, paste0(NAME, "_rf_fit.rds")))

  # temporary debug
  print(head(fit_glm$pred)); print(fit_glm$bestTune)


  # Pull pooled CV predictions and map rowIndex -> sample
  pred_glm <- best_cv_preds(fit_glm) |> mutate(model="glmnet")
  pred_rf  <- best_cv_preds(fit_rf)  |> mutate(model="rf")
  idx_map  <- tibble(rowIndex = seq_len(nrow(df_keep)), sample = df_keep$sample)

  pred_glm <- pred_glm |> left_join(idx_map, by="rowIndex")
  pred_rf  <- pred_rf  |> left_join(idx_map, by="rowIndex")

  write.csv(pred_glm, file.path(outdir, paste0(NAME, "_glmnet_oof_preds.csv")), row.names=FALSE)
  write.csv(pred_rf,  file.path(outdir, paste0(NAME, "_rf_oof_preds.csv")), row.names=FALSE)

  # Evaluate (AUC+CI, threshold, confusion)
  eval_from_preds(pred_glm, name=paste0(NAME, "_glmnet"), outdir=outdir)
  evrf <- eval_from_preds(pred_rf,  name=paste0(NAME, "_rf"),     outdir=outdir)

  # RF stability using best tune on full data
  rf_stability(X, y,
               bestTune = fit_rf$bestTune,
               runs = 50, seed0 = 2024,
               outdir = outdir,
               label = paste0(NAME, "_RF"))

  invisible(list(pred_glm=pred_glm, pred_rf=pred_rf))
}

# Run for all three
all_preds <- lapply(names(ACT_TABLES), eval_table)
names(all_preds) <- names(ACT_TABLES)

# ---------- Simple ensemble of GSVA+PROGENy (RF only) ----------
if (all(c("GSVA","PROGENy") %in% names(all_preds))) {
  ens_dir <- file.path(OUTDIR, "ENSEMBLE"); dir.create(ens_dir, showWarnings=FALSE)

  gs <- all_preds$GSVA$pred_rf |> select(sample, obs, prob) |> rename(prob_gsva=prob)
  pr <- all_preds$PROGENy$pred_rf |> select(sample, obs, prob) |> rename(prob_progeny=prob)

  ens <- full_join(gs, pr, by=c("sample","obs")) |>
    mutate(prob_mean = rowMeans(cbind(prob_gsva, prob_progeny), na.rm=TRUE))

  roc_ens <- pROC::roc(response = ens$obs,
                       predictor = ens$prob_mean,
                       levels = c("control","IOP"),
                       direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_ens))
  ci_auc  <- as.numeric(pROC::ci.auc(roc_ens))

  png(file.path(ens_dir, "Ensemble_RF_ROC.png"), width=1200, height=900, res=160)
  plot(roc_ens, main=paste0("RF Ensemble (GSVA+PROGENy) AUC=",
                            round(auc_val,3), " [", round(ci_auc[1],3), "–", round(ci_auc[3],3), "]"))
  abline(a=0,b=1,lty=2,col="grey60"); dev.off()

  write.csv(ens, file.path(ens_dir, "ensemble_oof_probs.csv"), row.names=FALSE)
  write.csv(tibble(model="RF_ensemble_GSVA_PROGENy",
                   auc=round(auc_val,3),
                   auc_ci_low=round(ci_auc[1],3),
                   auc_ci_high=round(ci_auc[3],3)),
            file.path(ens_dir, "ensemble_summary.csv"), row.names=FALSE)
}
