#!/usr/bin/env Rscript
#BiocManager::install("tidyverse")
#BiocManager::install("glmnet")
#BiocManager::install("ranger")
#BiocManager::install("pROC")

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(ranger)
  library(pROC)
})

meta_file <- "data/metadata/samples.txt"
tables <- list(
  DoRothEA = "results/pathways/DoRothEA_TFactivity.tsv",
  GSVA     = "results/pathways/GSVA_hallmark_scores.tsv",
  PROGENy  = "results/pathways/PROGENy_scores.tsv"
)
outdir <- "results/ml"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_meta <- function(f){
  m <- read.delim(f, header = TRUE, check.names = FALSE) %>% as_tibble()
  nms <- tolower(names(m))
  if(!all(c("sample","group") %in% nms)){
    if(ncol(m) >= 2) names(m)[1:2] <- c("sample","group")
    else stop("samples.txt must have at least 2 columns: sample, group")
  } else names(m) <- nms
  m %>% mutate(group = factor(group),
               y = factor(ifelse(group == levels(group)[1], "control", "IOP")))
}

read_activity <- function(path, samples){
  mat <- read.delim(path, header = TRUE, check.names = FALSE)
  if(!any(names(mat) == "feature") && !any(names(mat) == "gene")){
    if(all(!grepl("\\s", mat[[1]]))) {
      rn <- mat[[1]]; mat <- mat[,-1, drop=FALSE]; rownames(mat) <- rn
    }
  }
  cn <- colnames(mat); rn <- rownames(mat); smp <- samples$sample
  if(!is.null(cn) && length(intersect(cn, smp)) >= length(intersect(rn, smp))){
    X <- as_tibble(t(as.matrix(mat)), rownames = "sample")
  } else {
    X <- as_tibble(mat, rownames = "sample")
  }
  X <- X %>% filter(sample %in% smp)
  X <- samples %>% select(sample, y) %>% left_join(X, by="sample")
  # make numeric matrix
  feats <- X %>% select(-sample, -y) %>% mutate(across(everything(), as.numeric))
  list(X = X, y = X$y, feats = as.matrix(feats), featnames = colnames(feats))
}

stdz <- function(M){
  mu <- apply(M, 2, mean, na.rm = TRUE)
  sdv <- apply(M, 2, sd,   na.rm = TRUE); sdv[sdv==0 | !is.finite(sdv)] <- 1
  sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
}

make_folds <- function(y, k=5, seed=1){
  set.seed(seed)
  idx <- split(seq_along(y), y)
  folds <- vector("list", k)
  for(i in seq_len(k)){
    folds[[i]] <- unlist(lapply(idx, function(v) sample(v, ceiling(length(v)/k))))
    idx <- Map(setdiff, idx, list(folds[[i]]))
  }
  folds
}

plot_roc <- function(obs, prob, title){
  roc_obj <- pROC::roc(obs, prob, levels=c("control","IOP"), direction = "<")
  auc <- as.numeric(pROC::auc(roc_obj))
  p <- ggplot(data.frame(
      tpr = rev(roc_obj$sensitivities),
      fpr = rev(1 - roc_obj$specificities)
    ), aes(fpr, tpr)) +
    geom_line(color="#1b9e77", linewidth=1) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
    coord_equal() + theme_bw(base_size = 12) +
    labs(x="False positive rate", y="True positive rate",
         title=sprintf("%s (AUC = %.3f)", title, auc))
  list(plot = p, auc = auc)
}

for(tag in names(tables)){
  message("=== ", tag, " ===")
  meta <- read_meta(meta_file)
  dl <- read_activity(tables[[tag]], meta)
  X <- stdz(dl$feats)
  y <- dl$y
  featnames <- dl$featnames

  folds <- make_folds(y, k=5, seed=1)

  # ----------------- GLMNET (elastic net)
  oof_prob_glm <- rep(NA_real_, length(y))
  for(i in seq_along(folds)){
    te <- folds[[i]]; tr <- setdiff(seq_along(y), te)
    cv <- cv.glmnet(X[tr,], y[tr], family="binomial", alpha=0.5, type.measure="auc")
    mdl <- glmnet(X[tr,], y[tr], family="binomial", alpha=0.5, lambda=cv$lambda.1se)
    oof_prob_glm[te] <- drop(predict(mdl, X[te,], type="response"))
  }
  roc_glm <- plot_roc(y, oof_prob_glm, paste0(tag, " – glmnet"))

  ggsave(file.path(outdir, paste0(tag, "_glmnet_ROC.png")),
         roc_glm$plot, width=5, height=5, dpi=300)

  # final fit for feature importance
  cv_final <- cv.glmnet(X, y, family="binomial", alpha=0.5, type.measure="auc")
  mdl_final <- glmnet(X, y, family="binomial", alpha=0.5, lambda=cv_final$lambda.1se)
  coefs <- as.matrix(coef(mdl_final))[,1]
  imp_glm <- tibble(feature = rownames(as.matrix(coef(mdl_final))),
                    coef = as.numeric(coefs)) %>%
    filter(feature != "(Intercept)") %>%
    arrange(desc(abs(coef))) %>% slice_head(n=30)
  write.csv(imp_glm, file.path(outdir, paste0(tag, "_glmnet_feature_importance.csv")), row.names = FALSE)

  p_imp_glm <- imp_glm %>%
    mutate(feature = fct_reorder(feature, abs(coef))) %>%
    ggplot(aes(feature, coef, fill = coef > 0)) +
    geom_col() + coord_flip() + theme_bw(base_size = 11) +
    scale_fill_manual(values=c("#d95f02","#1b9e77")) +
    labs(title=paste0(tag, " – glmnet: top coefficients"), x=NULL, y="Coefficient")
  ggsave(file.path(outdir, paste0(tag, "_glmnet_feature_importance.png")),
         p_imp_glm, width=6, height=7, dpi=300)

  # ----------------- Random Forest (ranger)
  oof_prob_rf <- rep(NA_real_, length(y))
  for(i in seq_along(folds)){
    te <- folds[[i]]; tr <- setdiff(seq_along(y), te)
    rf <- ranger::ranger(
      y ~ ., data = data.frame(y=y[tr], X[tr,]),
      probability = TRUE, num.trees = 1000, mtry = max(1, floor(sqrt(ncol(X)))),
      importance = "impurity", seed = 1
    )
    phat <- predict(rf, data = data.frame(X[te,]))$predictions[, "IOP"]
    oof_prob_rf[te] <- phat
  }
  roc_rf <- plot_roc(y, oof_prob_rf, paste0(tag, " – random forest"))
  ggsave(file.path(outdir, paste0(tag, "_rf_ROC.png")),
         roc_rf$plot, width=5, height=5, dpi=300)

  # final RF fit for importance
  rf_final <- ranger::ranger(
    y ~ ., data = data.frame(y=y, X),
    probability = TRUE, num.trees = 1000,
    importance = "impurity", seed = 1
  )
  imp <- sort(rf_final$variable.importance, decreasing = TRUE)
  imp_rf <- tibble(feature = names(imp), importance = as.numeric(imp)) %>%
    slice_head(n=30)
  write.csv(imp_rf, file.path(outdir, paste0(tag, "_rf_feature_importance.csv")), row.names = FALSE)

  p_imp_rf <- imp_rf %>%
    mutate(feature = fct_reorder(feature, importance)) %>%
    ggplot(aes(feature, importance)) +
    geom_col(fill="#1b9e77") + coord_flip() + theme_bw(base_size = 11) +
    labs(title=paste0(tag, " – RF: top Gini importance"), x=NULL, y="Importance")
  ggsave(file.path(outdir, paste0(tag, "_rf_feature_importance.png")),
         p_imp_rf, width=6, height=7, dpi=300)
}

message("ML outputs written to: ", outdir)
