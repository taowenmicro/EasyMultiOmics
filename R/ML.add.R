# EasyMultiOmics biomarker ML v2
# Shared helpers + SVM/GLM/LASSO/XGBoost/Decision Tree/Bagging/Naive Bayes/NNet/ROC
# loadingPCA.metm2, LDA.metm2 and randomforest.metm2 are assumed to be defined separately.

.em_ml_prepare <- function(
    ps, group = "Group", rank = "Genus",
    prevalence = 0.30, min_mean = 0.01,
    variable_top = 20,
    transform = c("hellinger", "relative", "clr", "log1p", "none"),
    pseudocount = 0.5
) {
  transform <- match.arg(transform)
  ps <- phyloseq::filter_taxa(ps, function(x) sum(x) > 0, TRUE)

  if (!is.null(rank)) {
    tt <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
    if (is.null(tt)) stop("phyloseq对象没有tax_table，不能按 ", rank, " 聚合")
    tx <- as.data.frame(tt, stringsAsFactors = FALSE, check.names = FALSE)
    if (!rank %in% colnames(tx)) stop("tax_table中不存在分类层级: ", rank)
    z <- trimws(as.character(tx[[rank]]))
    ok <- !is.na(z) & nzchar(z) & !tolower(z) %in% c("unassigned", "unknown", "unclassified", "na")
    ps <- phyloseq::prune_taxa(rownames(tx)[ok], ps)
    if (phyloseq::ntaxa(ps) < 2) stop(rank, " 层级有效特征少于2个")
    ps <- phyloseq::tax_glom(ps, taxrank = rank, NArm = TRUE)
  }

  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)

  meta <- data.frame(phyloseq::sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE)
  if (!group %in% colnames(meta)) stop("sample_data中不存在分组列: ", group)
  meta <- meta[rownames(otu), , drop = FALSE]

  grp0 <- droplevels(factor(as.character(meta[[group]])))
  if (anyNA(grp0)) stop("分组变量存在NA")
  if (nlevels(grp0) != 2) stop("当前机器学习函数用于二元比较；检测到 ", nlevels(grp0), " 个分组")

  original_levels <- levels(grp0)
  safe_levels <- make.unique(make.names(original_levels))
  level_map <- setNames(safe_levels, original_levels)
  reverse_map <- setNames(original_levels, safe_levels)
  y <- factor(unname(level_map[as.character(grp0)]), levels = safe_levels)
  names(y) <- rownames(otu)

  lib <- rowSums(otu)
  if (any(lib <= 0)) stop("存在总reads为0的样本")
  rel <- sweep(otu, 1, lib, "/")
  rel[!is.finite(rel)] <- 0

  prev <- colMeans(otu > 0)
  mean_pct <- colMeans(rel) * 100
  eligible <- prev >= prevalence & mean_pct >= min_mean & apply(otu, 2, stats::sd, na.rm = TRUE) > 0
  if (!any(eligible)) stop("prevalence/min_mean过滤后无有效特征")

  X0 <- switch(
    transform,
    hellinger = sqrt(rel),
    relative = rel,
    clr = {
      z <- otu + pseudocount
      lz <- log(z)
      sweep(lz, 1, rowMeans(lz), "-")
    },
    log1p = log1p(rel * 1e6),
    none = otu
  )

  Xe <- X0[, eligible, drop = FALSE]
  variability <- apply(Xe, 2, stats::mad, na.rm = TRUE)
  if (!any(is.finite(variability) & variability > 0))
    variability <- apply(Xe, 2, stats::IQR, na.rm = TRUE)
  if (!any(is.finite(variability) & variability > 0))
    variability <- apply(Xe, 2, stats::var, na.rm = TRUE)

  selected <- colnames(Xe)
  if (!is.null(variable_top) && variable_top > 0 && length(selected) > variable_top) {
    ord <- order(variability, decreasing = TRUE, na.last = NA)
    selected <- names(variability)[head(ord, variable_top)]
  }

  X <- X0[, selected, drop = FALSE]
  sdv <- apply(X, 2, stats::sd, na.rm = TRUE)
  X <- X[, is.finite(sdv) & sdv > 1e-12, drop = FALSE]
  if (ncol(X) < 2) stop("无监督过滤后有效特征少于2个")

  tt <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
  tax <- if (is.null(tt)) data.frame(FeatureID = colnames(otu), stringsAsFactors = FALSE) else {
    q <- as.data.frame(tt, stringsAsFactors = FALSE, check.names = FALSE)
    q$FeatureID <- rownames(q)
    q
  }

  filter_table <- data.frame(
    FeatureID = colnames(otu),
    Prevalence = prev * 100,
    MeanAbundance = mean_pct,
    Variability = NA_real_,
    Eligible = eligible,
    Selected = colnames(otu) %in% colnames(X),
    stringsAsFactors = FALSE
  )
  ii <- match(names(variability), filter_table$FeatureID)
  filter_table$Variability[ii] <- variability

  list(
    ps = ps, X = as.data.frame(X, check.names = FALSE), y = y,
    y_original = grp0, rel = rel, otu = otu, meta = meta, tax = tax,
    filter_table = filter_table, level_map = level_map, reverse_map = reverse_map,
    group = group, rank = rank, transform = transform
  )
}

.em_model_colors <- function(groups, group_colors = NULL) {
  base <- setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
  if (is.null(group_colors)) return(base)
  z <- as.character(group_colors)
  if (is.null(names(z))) names(z) <- groups[seq_len(min(length(z), length(groups)))]
  miss <- setdiff(groups, names(z))
  c(z, base[miss])
}

.em_eval_pred <- function(pred, positive, negative, reverse_map) {
  pred$obs <- factor(pred$obs, levels = c(positive, negative))
  pred$pred <- factor(pred$pred, levels = c(positive, negative))
  cm <- caret::confusionMatrix(pred$pred, pred$obs, positive = positive)
  roc_obj <- tryCatch(
    pROC::roc(pred$obs, pred[[positive]], levels = c(negative, positive), direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  auc <- if (is.null(roc_obj)) NA_real_ else as.numeric(pROC::auc(roc_obj))
  by <- cm$byClass
  metric <- data.frame(
    Accuracy = unname(cm$overall["Accuracy"]),
    BalancedAccuracy = unname(by["Balanced Accuracy"]),
    Sensitivity = unname(by["Sensitivity"]),
    Specificity = unname(by["Specificity"]),
    AUC = auc
  )
  pred$Observed <- unname(reverse_map[as.character(pred$obs)])
  pred$Predicted <- unname(reverse_map[as.character(pred$pred)])
  list(metric = metric, confusion = as.matrix(cm$table), roc = roc_obj, pred = pred)
}

.em_importance_table <- function(fit, prep, lab = NULL, fill = "Phylum", top = 20) {
  vi <- tryCatch(caret::varImp(fit, scale = FALSE)$importance, error = function(e) NULL)
  if (is.null(vi) || !nrow(vi)) return(data.frame())
  vi <- as.data.frame(vi, check.names = FALSE)
  num <- vapply(vi, is.numeric, logical(1))
  vi$Importance <- if (any(num)) rowMeans(abs(as.matrix(vi[, num, drop = FALSE])), na.rm = TRUE) else NA_real_
  vi$FeatureID <- rownames(vi)
  vi <- dplyr::left_join(vi, prep$tax, by = "FeatureID")
  vi <- dplyr::left_join(vi, prep$filter_table[, c("FeatureID", "Prevalence", "MeanAbundance", "Variability")], by = "FeatureID")
  if (is.null(lab)) lab <- if (!is.null(prep$rank)) prep$rank else "FeatureID"
  vi$Label <- if (lab %in% colnames(vi)) as.character(vi[[lab]]) else vi$FeatureID
  bad <- is.na(vi$Label) | !nzchar(vi$Label) | tolower(vi$Label) %in% c("unassigned", "unknown", "unclassified", "na")
  vi$Label[bad] <- vi$FeatureID[bad]
  vi$Fill <- if (!is.null(fill) && fill %in% colnames(vi)) as.character(vi[[fill]]) else "Taxon"
  vi$Fill[is.na(vi$Fill) | !nzchar(vi$Fill)] <- "Unassigned"
  vi <- vi[order(-vi$Importance), , drop = FALSE]
  vi$Top <- seq_len(nrow(vi)) <= min(top, nrow(vi))
  vi
}

.em_caret_metm2 <- function(
    ps, model_name, method, group = "Group", rank = "Genus",
    k = 3, seed = 1010, top = 20,
    prevalence = 0.30, min_mean = 0.01, variable_top = 20,
    transform = "hellinger", lab = NULL, fill = "Phylum",
    group_colors = NULL, tuneGrid = NULL, tuneLength = NULL,
    preProcess = NULL, extra_args = list()
) {
  prep <- .em_ml_prepare(ps, group, rank, prevalence, min_mean, variable_top, transform)
  X <- prep$X; y <- prep$y
  k <- max(2L, min(as.integer(k), as.integer(min(table(y)))))
  if (min(table(y)) < 2) stop("每组至少需要2个样本")

  ctrl <- caret::trainControl(
    method = "cv", number = k, classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final", returnResamp = "final",
    allowParallel = FALSE
  )

  args <- list(
    x = X, y = y, method = method,
    metric = "ROC", trControl = ctrl
  )
  if (!is.null(tuneGrid)) args$tuneGrid <- tuneGrid
  if (is.null(tuneGrid) && !is.null(tuneLength)) args$tuneLength <- tuneLength
  if (!is.null(preProcess)) args$preProcess <- preProcess
  args <- c(args, extra_args)

  set.seed(seed)
  fit <- do.call(caret::train, args)

  pred <- fit$pred
  if (nrow(fit$bestTune) && nrow(pred)) {
    for (nm in intersect(names(fit$bestTune), names(pred))) {
      bt <- fit$bestTune[[nm]][1]
      if (is.numeric(pred[[nm]])) pred <- pred[abs(pred[[nm]] - bt) < 1e-12, , drop = FALSE]
      else pred <- pred[as.character(pred[[nm]]) == as.character(bt), , drop = FALSE]
    }
  }

  positive <- levels(y)[1]; negative <- levels(y)[2]
  ev <- .em_eval_pred(pred, positive, negative, prep$reverse_map)
  metric <- cbind(
    data.frame(Model = model_name, Samples = nrow(X), Features = ncol(X), CV_folds = k),
    ev$metric
  )

  imp <- .em_importance_table(fit, prep, lab, fill, top)
  top_imp <- head(imp, min(top, nrow(imp)))
  cols <- .em_model_colors(levels(prep$y_original), group_colors)

  p_imp <- NULL
  if (nrow(top_imp)) {
    top_imp$LabelPlot <- factor(make.unique(top_imp$Label), levels = rev(make.unique(top_imp$Label)))
    fill_levels <- unique(top_imp$Fill)
    fill_cols <- setNames(grDevices::hcl.colors(length(fill_levels), "Dynamic"), fill_levels)
    p_imp <- ggplot2::ggplot(top_imp, ggplot2::aes(Importance, LabelPlot, color = Fill)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = Importance, yend = LabelPlot), linewidth = 1.1) +
      ggplot2::geom_point(size = 3.5) +
      ggplot2::scale_color_manual(values = fill_cols) +
      ggplot2::labs(x = "Variable importance", y = NULL, color = fill) +
      ggplot2::theme_classic()
  }

  p_roc <- NULL; roc_data <- NULL
  if (!is.null(ev$roc)) {
    roc_data <- data.frame(
      FPR = 1 - ev$roc$specificities,
      TPR = ev$roc$sensitivities
    )
    p_roc <- ggplot2::ggplot(roc_data, ggplot2::aes(FPR, TPR)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = "1 - Specificity", y = "Sensitivity",
        title = paste0(model_name, " CV ROC, AUC = ", round(metric$AUC, 3))
      ) + ggplot2::theme_classic()
  }

  cm_long <- as.data.frame(as.table(ev$confusion))
  colnames(cm_long) <- c("PredictedSafe", "ActualSafe", "N")
  cm_long$Predicted <- unname(prep$reverse_map[as.character(cm_long$PredictedSafe)])
  cm_long$Actual <- unname(prep$reverse_map[as.character(cm_long$ActualSafe)])
  p_conf <- ggplot2::ggplot(cm_long, ggplot2::aes(Predicted, Actual, fill = N)) +
    ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = N), size = 4) +
    ggplot2::scale_fill_viridis_c() + ggplot2::theme_classic()

  prob_col <- positive
  p_prob <- ggplot2::ggplot(ev$pred, ggplot2::aes(Observed, .data[[prob_col]], color = Observed)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55) +
    ggplot2::geom_jitter(width = 0.08, size = 2) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(x = NULL, y = paste0("CV probability of ", prep$reverse_map[[positive]])) +
    ggplot2::theme_classic() + ggplot2::theme(legend.position = "none")

  list(
    metric, imp,
    importance_plot = p_imp, roc_plot = p_roc,
    confusion_plot = p_conf, probability_plot = p_prob,
    model = fit, predictions = ev$pred,
    confusion = ev$confusion, roc_data = roc_data,
    filter_table = prep$filter_table,
    parameters = data.frame(
      Parameter = c("model", "group", "rank", "k", "seed", "prevalence", "min_mean_percent", "variable_top", "transform"),
      Value = c(model_name, group, ifelse(is.null(rank), "NULL", rank), k, seed, prevalence, min_mean, variable_top, transform),
      stringsAsFactors = FALSE
    )
  )
}

svm_metm2 <- function(ps, group = "Group", rank = "Genus", k = 3, seed = 1010, top = 20,
                      prevalence = 0.30, min_mean = 0.01, variable_top = 10,
                      transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL) {
  .em_caret_metm2(
    ps, "SVM", "svmLinear2", group, rank, k, seed, top,
    prevalence, min_mean, variable_top, transform, lab, fill, group_colors,
    tuneGrid = data.frame(cost = 2^(-2:2)), preProcess = c("center", "scale")
  )
}

glm.metm2 <- function(ps, group = "Group", rank = "Genus", k = 3, seed = 1010, top = 20,
                      prevalence = 0.30, min_mean = 0.01, variable_top = 3,
                      transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL) {
  .em_caret_metm2(
    ps, "GLM", "glm", group, rank, k, seed, top,
    prevalence, min_mean, variable_top, transform, lab, fill, group_colors,
    extra_args = list(family = stats::binomial())
  )
}

lasso.metm2 <- function(ps, group = "Group", rank = "Genus", k = 3, seed = 1010, top = 20,
                        prevalence = 0.30, min_mean = 0.01, variable_top = 20,
                        transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL) {
  grid <- expand.grid(alpha = 1, lambda = 10^seq(-4, 0, length.out = 25))
  .em_caret_metm2(
    ps, "LASSO", "glmnet", group, rank, k, seed, top,
    prevalence, min_mean, variable_top, transform, lab, fill, group_colors,
    tuneGrid = grid, preProcess = c("center", "scale")
  )
}

xgboost.metm2 <- function(ps, group = "Group", rank = "Genus", k = 3, seed = 1010, top = 20,
                          prevalence = 0.30, min_mean = 0.01, variable_top = 20,
                          transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL) {
  grid <- expand.grid(
    nrounds = c(25, 50), max_depth = c(1, 2), eta = c(0.05, 0.1),
    gamma = 0, colsample_bytree = 1, min_child_weight = 1, subsample = 0.8
  )
  .em_caret_metm2(
    ps, "XGBoost", "xgbTree", group, rank, k, seed, top,
    prevalence, min_mean, variable_top, transform, lab, fill, group_colors,
    tuneGrid = grid, extra_args = list(verbose = 0, nthread = 1)
  )
}

decisiontree.metm2 <- function(ps, group = "Group", rank = "Genus", k = 3, seed = 6358, top = 20,
                               prevalence = 0.30, min_mean = 0.01, variable_top = 20,
                               transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL) {
  .em_caret_metm2(
    ps, "DecisionTree", "rpart", group, rank, k, seed, top,
    prevalence, min_mean, variable_top, transform, lab, fill, group_colors,
    tuneGrid = data.frame(cp = c(0, 0.01, 0.05, 0.1))
  )
}

bagging_metm2 <- function(ps, group = "Group", rank = "Genus", k = 3, seed = 1010, top = 20,
                          prevalence = 0.30, min_mean = 0.01, variable_top = 20,
                          transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL) {
  .em_caret_metm2(
    ps, "Bagging", "treebag", group, rank, k, seed, top,
    prevalence, min_mean, variable_top, transform, lab, fill, group_colors
  )
}

nnet.metm2 <- function(ps, group = "Group", rank = "Genus", k = 3, seed = 1010, top = 20,
                       prevalence = 0.30, min_mean = 0.01, variable_top = 10,
                       transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL) {
  grid <- expand.grid(size = c(1, 2, 3), decay = c(0, 0.01, 0.1))
  .em_caret_metm2(
    ps, "NNet", "nnet", group, rank, k, seed, top,
    prevalence, min_mean, variable_top, transform, lab, fill, group_colors,
    tuneGrid = grid, preProcess = c("center", "scale"),
    extra_args = list(trace = FALSE, MaxNWts = 50000, maxit = 500)
  )
}

naivebayes.metm2 <- function(
    ps, group = "Group", rank = "Genus", k = 3, seed = 1010, top = 20,
    prevalence = 0.30, min_mean = 0.01, variable_top = 10,
    transform = "hellinger", lab = NULL, fill = "Phylum", group_colors = NULL,
    laplace = 1
) {
  if (!requireNamespace("e1071", quietly = TRUE)) stop("naivebayes.metm2需要e1071")
  prep <- .em_ml_prepare(ps, group, rank, prevalence, min_mean, variable_top, transform)
  X <- prep$X; y <- prep$y
  k <- max(2L, min(as.integer(k), as.integer(min(table(y)))))
  set.seed(seed)
  folds <- caret::createFolds(y, k = k, returnTrain = FALSE)
  positive <- levels(y)[1]; negative <- levels(y)[2]

  pred_list <- lapply(seq_along(folds), function(i) {
    te <- folds[[i]]; tr <- setdiff(seq_len(nrow(X)), te)
    fit <- e1071::naiveBayes(x = X[tr, , drop = FALSE], y = y[tr], laplace = laplace)
    pr <- predict(fit, X[te, , drop = FALSE], type = "raw")
    cl <- predict(fit, X[te, , drop = FALSE], type = "class")
    data.frame(rowIndex = te, obs = y[te], pred = cl, Prob = pr[, positive], Resample = paste0("Fold", i))
  })
  pred <- dplyr::bind_rows(pred_list)
  pred[[positive]] <- pred$Prob
  ev <- .em_eval_pred(pred, positive, negative, prep$reverse_map)

  full_model <- e1071::naiveBayes(x = X, y = y, laplace = laplace)
  g1 <- levels(y)[1]; g2 <- levels(y)[2]
  sep <- vapply(seq_len(ncol(X)), function(j) {
    a <- X[y == g1, j]; b <- X[y == g2, j]
    abs(mean(a) - mean(b)) / sqrt((stats::var(a) + stats::var(b)) / 2 + 1e-12)
  }, numeric(1))
  imp <- data.frame(FeatureID = colnames(X), Importance = sep, stringsAsFactors = FALSE)
  imp <- dplyr::left_join(imp, prep$tax, by = "FeatureID")
  imp <- dplyr::left_join(imp, prep$filter_table[, c("FeatureID", "Prevalence", "MeanAbundance", "Variability")], by = "FeatureID")
  if (is.null(lab)) lab <- if (!is.null(rank)) rank else "FeatureID"
  imp$Label <- if (lab %in% colnames(imp)) as.character(imp[[lab]]) else imp$FeatureID
  bad <- is.na(imp$Label) | !nzchar(imp$Label); imp$Label[bad] <- imp$FeatureID[bad]
  imp$Fill <- if (!is.null(fill) && fill %in% colnames(imp)) as.character(imp[[fill]]) else "Taxon"
  imp$Fill[is.na(imp$Fill) | !nzchar(imp$Fill)] <- "Unassigned"
  imp <- imp[order(-imp$Importance), , drop = FALSE]

  metric <- cbind(data.frame(Model = "NaiveBayes", Samples = nrow(X), Features = ncol(X), CV_folds = k), ev$metric)
  top_imp <- head(imp, min(top, nrow(imp)))
  p_imp <- NULL
  if (nrow(top_imp)) {
    top_imp$LabelPlot <- factor(make.unique(top_imp$Label), levels = rev(make.unique(top_imp$Label)))
    p_imp <- ggplot2::ggplot(top_imp, ggplot2::aes(Importance, LabelPlot)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = Importance, yend = LabelPlot), linewidth = 1.1) +
      ggplot2::geom_point(size = 3.5) + ggplot2::labs(x = "NB class-separation score", y = NULL) + ggplot2::theme_classic()
  }
  p_roc <- NULL; roc_data <- NULL
  if (!is.null(ev$roc)) {
    roc_data <- data.frame(FPR = 1 - ev$roc$specificities, TPR = ev$roc$sensitivities)
    p_roc <- ggplot2::ggplot(roc_data, ggplot2::aes(FPR, TPR)) + ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + ggplot2::coord_equal() +
      ggplot2::labs(title = paste0("Naive Bayes CV ROC, AUC = ", round(metric$AUC, 3))) + ggplot2::theme_classic()
  }
  cm_long <- as.data.frame(as.table(ev$confusion)); colnames(cm_long) <- c("PredictedSafe", "ActualSafe", "N")
  cm_long$Predicted <- unname(prep$reverse_map[as.character(cm_long$PredictedSafe)])
  cm_long$Actual <- unname(prep$reverse_map[as.character(cm_long$ActualSafe)])
  p_conf <- ggplot2::ggplot(cm_long, ggplot2::aes(Predicted, Actual, fill = N)) + ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = N), size = 4) + ggplot2::scale_fill_viridis_c() + ggplot2::theme_classic()

  list(
    metric, imp, importance_plot = p_imp, roc_plot = p_roc, confusion_plot = p_conf,
    model = full_model, predictions = ev$pred, confusion = ev$confusion,
    roc_data = roc_data, filter_table = prep$filter_table,
    parameters = data.frame(Parameter = c("group", "rank", "k", "seed", "prevalence", "min_mean_percent", "variable_top", "transform", "laplace"),
                            Value = c(group, rank, k, seed, prevalence, min_mean, variable_top, transform, laplace), stringsAsFactors = FALSE)
  )
}

Roc.metm2 <- function(
    ps, group = "Group", rank = "Genus", top = 10,
    prevalence = 0.30, min_mean = 0.01,
    boot_n = 500, seed = 11, group_colors = NULL
) {
  prep <- .em_ml_prepare(ps, group, rank, prevalence, min_mean, variable_top = NULL, transform = "relative")
  rel <- prep$rel[, colnames(prep$X), drop = FALSE] * 100
  y <- prep$y_original
  lev <- levels(y)
  set.seed(seed)

  rows <- vector("list", ncol(rel)); curves <- vector("list", ncol(rel))
  for (j in seq_len(ncol(rel))) {
    id <- colnames(rel)[j]; x <- rel[, j]
    rr <- tryCatch(pROC::roc(y, x, levels = lev, direction = "auto", quiet = TRUE), error = function(e) NULL)
    if (is.null(rr)) next
    ci <- tryCatch(as.numeric(pROC::ci.auc(rr, method = "bootstrap", boot.n = boot_n, stratified = TRUE)), error = function(e) c(NA, NA, NA))
    best <- tryCatch(pROC::coords(rr, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE), error = function(e) NULL)
    rows[[j]] <- data.frame(
      FeatureID = id, AUC = as.numeric(pROC::auc(rr)), CI_low = ci[1], CI_mid = ci[2], CI_high = ci[3],
      Threshold = if (is.null(best)) NA_real_ else as.numeric(best$threshold),
      Sensitivity = if (is.null(best)) NA_real_ else as.numeric(best$sensitivity),
      Specificity = if (is.null(best)) NA_real_ else as.numeric(best$specificity),
      Mean_Group1 = mean(x[y == lev[1]]), Mean_Group2 = mean(x[y == lev[2]]), stringsAsFactors = FALSE
    )
    curves[[j]] <- data.frame(FeatureID = id, FPR = 1 - rr$specificities, TPR = rr$sensitivities)
  }
  tab <- dplyr::bind_rows(rows)
  curves <- dplyr::bind_rows(curves)
  tab <- dplyr::left_join(tab, prep$tax, by = "FeatureID")
  tab <- tab[order(-tab$AUC), , drop = FALSE]
  lab <- if (!is.null(rank) && rank %in% names(tab)) rank else "FeatureID"
  tab$Label <- as.character(tab[[lab]]); bad <- is.na(tab$Label) | !nzchar(tab$Label); tab$Label[bad] <- tab$FeatureID[bad]
  top_tab <- head(tab, min(top, nrow(tab)))
  top_ids <- top_tab$FeatureID
  curve_top <- curves[curves$FeatureID %in% top_ids, , drop = FALSE]
  curve_top$Label <- setNames(make.unique(top_tab$Label), top_tab$FeatureID)[curve_top$FeatureID]

  p_roc <- ggplot2::ggplot(curve_top, ggplot2::aes(FPR, TPR, color = Label)) +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::coord_equal() + ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", color = rank) + ggplot2::theme_classic()

  top_tab$LabelPlot <- factor(make.unique(top_tab$Label), levels = rev(make.unique(top_tab$Label)))
  p_auc <- ggplot2::ggplot(top_tab, ggplot2::aes(AUC, LabelPlot)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0.5, xend = AUC, yend = LabelPlot), linewidth = 1) +
    ggplot2::geom_point(size = 3.5) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2) + ggplot2::coord_cartesian(xlim = c(0.5, 1)) +
    ggplot2::labs(x = "AUC (bootstrap 95% CI)", y = NULL) + ggplot2::theme_classic()

  ab <- as.data.frame(rel[, top_ids, drop = FALSE], check.names = FALSE)
  ab$SampleID <- rownames(ab); ab$Group <- as.character(y)
  ab <- tidyr::pivot_longer(ab, cols = -c(SampleID, Group), names_to = "FeatureID", values_to = "Abundance")
  ab$Label <- setNames(make.unique(top_tab$Label), top_tab$FeatureID)[ab$FeatureID]
  cols <- .em_model_colors(lev, group_colors)
  p_ab <- ggplot2::ggplot(ab, ggplot2::aes(Group, Abundance, color = Group)) +
    ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.08, size = 1.8) +
    ggplot2::facet_wrap(~Label, scales = "free_y", ncol = 3) + ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(y = "Relative abundance (%)") + ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  list(
    p_roc, tab, curves,
    auc_plot = p_auc, abundance_plot = p_ab,
    top_results = top_tab, abundance_data = ab,
    parameters = data.frame(Parameter = c("group", "rank", "top", "prevalence", "min_mean_percent", "boot_n", "seed"),
                            Value = c(group, rank, top, prevalence, min_mean, boot_n, seed), stringsAsFactors = FALSE)
  )
}

# Compatibility wrapper: the new randomforest.metm2 already performs RF-CV.
rfcv.metm2 <- function(ps, group = "Group", rank = "Genus", nrfcvnum = 3, ...) {
  z <- randomforest.metm2(ps = ps, group = group, rank = rank, rfcv = TRUE, nrfcvnum = nrfcvnum, ...)
  list(z$rfcv_plot, z$rfcv_data, z$model_stats, randomforest_result = z)
}


# ============================================================================
# Pairwise PERMANOVA + workflow helpers (v3)
# ============================================================================

# Safe fallback helpers. If the user's pipeline already defines these, keep the
# existing definitions.
if (!exists("is_error_result", mode = "function")) {
  is_error_result <- function(x) {
    inherits(x, "em_error_result") || inherits(x, "try-error") ||
      inherits(x, "error") ||
      (is.list(x) && isTRUE(x$error))
  }
}

if (!exists("run_safe", mode = "function")) {
  run_safe <- function(name, expr) {
    tryCatch(
      expr,
      error = function(e) {
        message("[ERROR] ", name, ": ", conditionMessage(e))
        structure(
          list(error = TRUE, name = name, message = conditionMessage(e)),
          class = "em_error_result"
        )
      }
    )
  }
}

# Replace an existing worksheet case-insensitively, then write the new table.
write_sheet_replace <- function(wb, sheet, dat) {
  if (is.null(dat)) return(invisible(NULL))

  dat <- tryCatch(
    {
      if (exists("as_df_safe", mode = "function")) as_df_safe(dat) else as.data.frame(dat)
    },
    error = function(e) as.data.frame(dat)
  )

  if (!nrow(dat) && !ncol(dat)) return(invisible(NULL))

  sheet <- substr(as.character(sheet), 1, 31)
  old <- names(wb)
  hit <- which(tolower(old) == tolower(sheet))
  if (length(hit)) openxlsx::removeWorksheet(wb, old[hit[1]])

  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, dat)
  invisible(NULL)
}

# Pairwise PERMANOVA intended to decide which group pairs proceed to supervised
# pairwise ML. Bray-Curtis defaults to relative abundance.
pairwise_permanova_ml <- function(
    ps,
    group = "Group",
    dist = "bray",
    transform = c("relative", "none", "hellinger"),
    permutations = 999,
    seed = 11,
    adjust.method = "BH",
    selection = c("padj", "raw", "exploratory_3x3"),
    p_cutoff = 0.05,
    padj_cutoff = 0.05,
    min_R2 = 0,
    check_dispersion = TRUE,
    require_homogeneous_dispersion = FALSE,
    dispersion_alpha = 0.05
) {
  transform <- match.arg(transform)
  selection <- match.arg(selection)

  meta0 <- data.frame(
    phyloseq::sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!group %in% names(meta0)) stop("sample_data中不存在分组列: ", group)

  groups <- sort(unique(as.character(meta0[[group]])))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) < 2) stop("至少需要2个分组")

  pairs <- utils::combn(groups, 2, simplify = FALSE)

  out <- lapply(seq_along(pairs), function(ii) {
    g <- pairs[[ii]]

    ids <- rownames(meta0)[as.character(meta0[[group]]) %in% g]
    ps2 <- phyloseq::prune_samples(ids, ps)
    ps2 <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps2) > 0, ps2)

    meta2 <- data.frame(
      phyloseq::sample_data(ps2),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    meta2[[group]] <- droplevels(factor(as.character(meta2[[group]]), levels = g))

    n1 <- sum(meta2[[group]] == g[1], na.rm = TRUE)
    n2 <- sum(meta2[[group]] == g[2], na.rm = TRUE)

    if (n1 < 2 || n2 < 2) {
      return(data.frame(
        Group1 = g[1], Group2 = g[2], N1 = n1, N2 = n2,
        R2 = NA_real_, F = NA_real_, P = NA_real_,
        Dispersion_F = NA_real_, Dispersion_P = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    psd <- switch(
      transform,
      relative = phyloseq::transform_sample_counts(ps2, function(x) x / sum(x)),
      hellinger = phyloseq::transform_sample_counts(ps2, function(x) sqrt(x / sum(x))),
      none = ps2
    )

    d <- phyloseq::distance(psd, method = dist)
    labs <- attr(d, "Labels")
    meta2 <- meta2[labs, , drop = FALSE]

    set.seed(seed + ii - 1L)
    ad <- vegan::adonis2(
      d ~ .group_tmp,
      data = data.frame(.group_tmp = meta2[[group]], row.names = rownames(meta2)),
      permutations = permutations
    )

    disp_F <- disp_P <- NA_real_
    if (check_dispersion) {
      bd <- tryCatch(
        vegan::betadisper(d, meta2[[group]]),
        error = function(e) NULL
      )
      if (!is.null(bd)) {
        set.seed(seed + 10000L + ii - 1L)
        dp <- tryCatch(
          vegan::permutest(bd, permutations = permutations),
          error = function(e) NULL
        )
        if (!is.null(dp) && !is.null(dp$tab)) {
          disp_F <- suppressWarnings(as.numeric(dp$tab[1, "F"]))
          disp_P <- suppressWarnings(as.numeric(dp$tab[1, "Pr(>F)"]))
        }
      }
    }

    data.frame(
      Group1 = g[1],
      Group2 = g[2],
      N1 = n1,
      N2 = n2,
      R2 = as.numeric(ad$R2[1]),
      F = as.numeric(ad$F[1]),
      P = as.numeric(ad$`Pr(>F)`[1]),
      Dispersion_F = disp_F,
      Dispersion_P = disp_P,
      stringsAsFactors = FALSE
    )
  })

  tab <- dplyr::bind_rows(out)
  tab$Padj <- stats::p.adjust(tab$P, method = adjust.method)

  # Approximate lower bound for equal-size 3-vs-3 label permutations.
  # For n1 == n2, swapping group labels gives the same pseudo-F, so the
  # smallest attainable exact permutation P is at least 2 / choose(n1+n2,n1).
  tab$Approx_min_P <- ifelse(
    tab$N1 == tab$N2 & is.finite(tab$N1),
    2 / choose(tab$N1 + tab$N2, tab$N1),
    1 / choose(tab$N1 + tab$N2, tab$N1)
  )

  pass_disp <- rep(TRUE, nrow(tab))
  if (require_homogeneous_dispersion) {
    pass_disp <- is.na(tab$Dispersion_P) | tab$Dispersion_P > dispersion_alpha
  }

  tab$ML_selected <- switch(
    selection,
    padj = !is.na(tab$Padj) & tab$Padj <= padj_cutoff & tab$R2 >= min_R2 & pass_disp,
    raw = !is.na(tab$P) & tab$P <= p_cutoff & tab$R2 >= min_R2 & pass_disp,
    exploratory_3x3 = !is.na(tab$P) & tab$P <= max(0.10, p_cutoff) &
      tab$R2 >= min_R2 & pass_disp
  )

  tab$Selection_mode <- selection
  tab <- tab[order(tab$Padj, tab$P, -tab$R2), , drop = FALSE]
  rownames(tab) <- NULL

  if (selection == "padj" && any(tab$N1 == 3 & tab$N2 == 3, na.rm = TRUE) && padj_cutoff <= 0.05) {
    message(
      "注意：3 vs 3 的两组置换PERMANOVA中，因可交换标签数量极少，",
      "最小原始P通常约为0.10；因此 Padj <= 0.05 基本不可能达到。",
      "正式分析建议保留这一限制；若仅做探索，可显式使用 selection='exploratory_3x3'。"
    )
  }

  tab
}

# Save a standard pairwise ML result into its own model folder and pair workbook.
save_pair_ml_result <- function(
    res,
    tag,
    pair_path,
    pair_wb,
    pair_xlsx,
    width = 10,
    height = 8
) {
  if (is.null(res) || is_error_result(res)) return(invisible(NULL))

  model_path <- file.path(pair_path, tag)
  dir.create(model_path, recursive = TRUE, showWarnings = FALSE)

  if (length(res) >= 1 && is.data.frame(res[[1]]) && nrow(res[[1]]))
    write_sheet_replace(pair_wb, paste0(tag, "_metric"), res[[1]])

  if (length(res) >= 2 && is.data.frame(res[[2]]) && nrow(res[[2]]))
    write_sheet_replace(pair_wb, paste0(tag, "_importance"), res[[2]])

  if (!is.null(res$filter_table))
    write_sheet_replace(pair_wb, paste0(tag, "_filter"), res$filter_table)
  if (!is.null(res$predictions))
    write_sheet_replace(pair_wb, paste0(tag, "_pred"), res$predictions)
  if (!is.null(res$parameters))
    write_sheet_replace(pair_wb, paste0(tag, "_param"), res$parameters)
  if (!is.null(res$model_stats))
    write_sheet_replace(pair_wb, paste0(tag, "_stats"), res$model_stats)
  if (!is.null(res$confusion))
    write_sheet_replace(pair_wb, paste0(tag, "_confusion"), res$confusion)
  if (!is.null(res$class_performance))
    write_sheet_replace(pair_wb, paste0(tag, "_performance"), res$class_performance)
  if (!is.null(res$rfcv_data))
    write_sheet_replace(pair_wb, paste0(tag, "_featureCV"), res$rfcv_data)
  if (!is.null(res$abundance_summary))
    write_sheet_replace(pair_wb, paste0(tag, "_abundance"), res$abundance_summary)

  plots <- list(
    importance = res$importance_plot,
    roc = res$roc_plot,
    confusion = res$confusion_plot,
    probability = res$probability_plot,
    variability = res$variability_plot,
    error = res$error_plot,
    mds = res$mds_plot,
    heatmap = res$abundance_heatmap,
    abundance = res$abundance_plot,
    combined = res$combined_plot,
    feature_cv = res$rfcv_plot
  )

  for (nm in names(plots)) {
    p <- plots[[nm]]
    if (is.null(p)) next
    if (!exists("save_plot2", mode = "function")) next
    tryCatch(
      save_plot2(
        p, model_path, paste0(tolower(tag), "_", nm),
        width = if (nm == "combined") 18 else width,
        height = if (nm == "abundance") 12 else height
      ),
      error = function(e) message("[Plot skipped] ", tag, "/", nm, ": ", conditionMessage(e))
    )
  }

  openxlsx::saveWorkbook(pair_wb, pair_xlsx, overwrite = TRUE)
  invisible(res)
}
