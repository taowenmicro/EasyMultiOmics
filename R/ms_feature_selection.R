#' MSFeatureSelection: Comprehensive Untargeted Metabolomics Feature Selection
#'
#' 针对非靶向代谢组（MS）数据的特征挑选流水线：缺失与批次/漂移校正、归一化、QC 稳定性过滤、
#' 差异分析（limma）、机器学习重要性（RF/GBM/Boruta）、相关网络、综合加权排名。
#'
#' _PACKAGE
#' @name MSFeatureSelection
#' @importFrom utils write.csv
#' @importFrom sva ComBat
#' @importFrom Hmisc rcorr
#' @importFrom caret train trainControl varImp
NULL

# ========================= 1. 主函数 =========================

#' Untargeted Metabolomics Feature Selection (MS)
#'
#' @param ms_mat 矩阵/数据框：样本×代谢物（intensity），行名=样本名，列名=代谢物ID
#' @param sample_meta 数据框：包含分组/批次/QC 标识，行名=样本名
#' @param group_var 分组变量名（在 sample_meta 中）
#' @param batch_var 批次变量名（可选，用于ComBat）
#' @param qc_flag_var QC 标识变量名（TRUE/FALSE 或 0/1），用于漂移校正/稳定性过滤
#' @param prevalence_threshold 特征保留的样本占比阈值（默认0.7）
#' @param impute_method 缺失值填补：'halfmin'|'knn'（默认'halfmin'）
#' @param normalize_method 强度归一化：'TIC'|'median'|'quantile'（默认'median'）
#' @param combat 是否执行 ComBat 批次校正（默认 TRUE）
#' @param qc_drift_correction 是否执行 QC-based 漂移校正（默认 TRUE）
#' @param qc_rsd_cut 基于 QC 的 RSD 过滤阈值（默认0.3）
#' @param diff_method 差异分析方法：'limma'|'wilcoxon'（默认 'limma'）
#' @param ml_method ML 方法：'rf'|'gbm'|'boruta'（默认 'rf'）
#' @param cor_method 相关：'spearman'|'pearson'（默认 'spearman'）
#' @param cor_threshold 边阈值（默认 0.6）
#' @param p_threshold 相关 p 值阈值（默认 0.05）
#' @param weights 综合权重 list(abundance, importance, differential, network)
#' @param nfolds ML 交叉验证折数（默认 5）
#' @param ntree RF 森林数（默认 500）
#'
#' @return list: final_ranking / full_results / diff_results / importance_results /
#'   network_results / abundance_results / processed_ms（处理后矩阵）/ parameters
#' @export
ms_feature_selection <- function(ms_mat,
                                 sample_meta,
                                 group_var,
                                 batch_var = NULL,
                                 qc_flag_var = NULL,
                                 prevalence_threshold = 0.7,
                                 impute_method = c("halfmin","knn"),
                                 normalize_method = c("median","TIC","quantile"),
                                 combat = TRUE,
                                 qc_drift_correction = TRUE,
                                 qc_rsd_cut = 0.30,
                                 diff_method = c("limma","wilcoxon"),
                                 ml_method = c("rf","gbm","boruta"),
                                 cor_method = c("spearman","pearson"),
                                 cor_threshold = 0.6,
                                 p_threshold = 0.05,
                                 weights = list(abundance = 0.3,
                                                importance = 0.3,
                                                differential = 0.2,
                                                network = 0.2),
                                 nfolds = 5,
                                 ntree = 500) {

  impute_method   <- match.arg(impute_method)
  normalize_method<- match.arg(normalize_method)
  diff_method     <- match.arg(diff_method)
  ml_method       <- match.arg(ml_method)
  cor_method      <- match.arg(cor_method)

  stopifnot(is.matrix(ms_mat) || is.data.frame(ms_mat))
  ms_mat <- as.matrix(ms_mat)

  # ------- 样本对齐 -------
  common_samples <- intersect(rownames(ms_mat), rownames(sample_meta))
  if (length(common_samples) < 3) stop("样本行名需与元数据行名匹配，并且数量足够。")
  ms_mat <- ms_mat[common_samples, , drop = FALSE]
  sample_meta <- sample_meta[common_samples, , drop = FALSE]

  # 分组变量
  if (!group_var %in% colnames(sample_meta)) {
    stop(paste0("分组变量 ", group_var, " 不在 sample_meta 中"))
  }
  sample_meta[[group_var]] <- as.factor(sample_meta[[group_var]])

  # QC 标记/批次
  qc_vec <- NULL
  if (!is.null(qc_flag_var) && qc_flag_var %in% colnames(sample_meta)) {
    qc_vec <- as.logical(sample_meta[[qc_flag_var]])
    qc_vec[is.na(qc_vec)] <- FALSE
  }
  batch_vec <- NULL
  if (!is.null(batch_var) && batch_var %in% colnames(sample_meta)) {
    batch_vec <- as.factor(sample_meta[[batch_var]])
  }

  cat("==== MS feature selection start ====\n")
  cat(sprintf("Samples: %d, Features: %d\n", nrow(ms_mat), ncol(ms_mat)))
  cat(sprintf("Groups: %s\n", paste(levels(sample_meta[[group_var]]), collapse=", ")))
  if (!is.null(batch_vec)) cat(sprintf("Batches: %s\n", paste(levels(batch_vec), collapse=", ")))
  if (!is.null(qc_vec)) cat(sprintf("QC samples: %d\n", sum(qc_vec)))

  # ------- 1. 预处理（缺失过滤/补齐、归一化、批次/漂移） -------
  prep <- ms_preprocess(ms_mat, sample_meta,
                        prevalence_threshold = prevalence_threshold,
                        impute_method = impute_method,
                        normalize_method = normalize_method,
                        combat = combat,
                        batch_vec = batch_vec,
                        qc_flag = qc_vec,
                        qc_drift_correction = qc_drift_correction,
                        qc_rsd_cut = qc_rsd_cut)

  X <- prep$ms_proc           # 处理后的矩阵（样本×代谢物）
  sample_meta <- prep$sample_meta

  # ------- 2. 差异分析（limma / wilcoxon） -------
  cat("\n[Step 2] Differential analysis ...\n")
  diff_res <- ms_differential(X, sample_meta, group_var, method = diff_method)

  # ------- 3. 机器学习重要性 -------
  cat("\n[Step 3] Machine-learning importance ...\n")
  imp_res <- ms_ml_importance(X, sample_meta, group_var,
                              method = ml_method, nfolds = nfolds, ntree = ntree)

  # ------- 4. 相关网络 -------
  cat("\n[Step 4] Correlation network ...\n")
  net_res <- ms_network(X, cor_method = cor_method,
                        cor_threshold = cor_threshold, p_threshold = p_threshold)

  # ------- 5. 丰度/流行率（用于偏向高丰度） -------
  cat("\n[Step 5] Abundance metrics ...\n")
  abd_res <- ms_abundance_metrics(X)

  # ------- 6. 综合排名 -------
  cat("\n[Step 6] Integrated weighted ranking ...\n")
  rank_res <- ms_weighted_ranking(diff_res, imp_res, net_res, abd_res, weights)

  cat("\nDone.\n")
  cat(sprintf("Significant features (p<0.05): %d\n",
              sum(diff_res$pvalue < 0.05, na.rm = TRUE)))
  cat(sprintf("Network edges: %d\n", sum(net_res$Degree > 0)))

  res <- list(
    final_ranking = rank_res$summary_results,
    full_results  = rank_res$full_results,
    diff_results  = diff_res,
    importance_results = imp_res,
    network_results = net_res,
    abundance_results = abd_res,
    processed_ms  = X,
    parameters = list(
      group_var = group_var,
      prevalence_threshold = prevalence_threshold,
      impute_method = impute_method,
      normalize_method = normalize_method,
      combat = combat,
      qc_drift_correction = qc_drift_correction,
      qc_rsd_cut = qc_rsd_cut,
      diff_method = diff_method,
      ml_method = ml_method,
      cor_method = cor_method,
      cor_threshold = cor_threshold,
      p_threshold = p_threshold,
      weights = weights
    )
  )
  class(res) <- c("ms_feature_selection","list")
  return(res)
}

# ========================= 2. 预处理流水线 =========================

#' @keywords internal
ms_preprocess <- function(ms_mat, sample_meta,
                          prevalence_threshold = 0.7,
                          impute_method = "halfmin",
                          normalize_method = "median",
                          combat = TRUE,
                          batch_vec = NULL,
                          qc_flag = NULL,
                          qc_drift_correction = TRUE,
                          qc_rsd_cut = 0.30) {

  X <- as.matrix(ms_mat)

  # (a) 预过滤：按流行率（非缺失率）
  keep <- colMeans(is.finite(X) & X > 0, na.rm = TRUE) >= prevalence_threshold
  X <- X[, keep, drop = FALSE]
  cat(sprintf("Prevalence filter kept %d/%d features.\n", ncol(X), ncol(ms_mat)))

  # (b) 缺失填补
  X <- ms_impute(X, method = impute_method)

  # (c) 归一化（样本层）
  X <- ms_normalize(X, method = normalize_method)

  # (d) QC 漂移校正（如有 QC）
  if (qc_drift_correction && !is.null(qc_flag) && any(qc_flag)) {
    X <- ms_qc_drift_correction(X, qc_flag, sample_meta)
    # QC 稳定性过滤（基于QC CV/RSD）
    X <- ms_qc_rsd_filter(X, qc_flag, rsd_cut = qc_rsd_cut)
  }

  # (e) 批次校正（ComBat）
  if (combat && !is.null(batch_vec)) {
    if (!requireNamespace("sva", quietly = TRUE))
      stop("需要 sva::ComBat，请先安装 sva 包。")
    # ComBat 期望特征×样本
    Xc <- sva::ComBat(t(X), batch = batch_vec, par.prior = TRUE, prior.plots = FALSE)
    X <- t(Xc)
    cat("ComBat batch correction done.\n")
  }

  # (f) log2 变换（可选，这里默认所有值>0时再做）
  if (all(X > 0)) {
    X <- log2(X + 1)
  }

  # (g) 特征零方差过滤
  keep2 <- apply(X, 2, stats::var, na.rm = TRUE) > 1e-12
  X <- X[, keep2, drop = FALSE]

  return(list(ms_proc = X, sample_meta = sample_meta))
}

# ---------- 缺失填补 ----------
#' @keywords internal
ms_impute <- function(X, method = "halfmin") {
  if (method == "halfmin") {
    x_min <- apply(X, 2, function(v) min(v[v>0], na.rm = TRUE))
    x_min[!is.finite(x_min)] <- 1e-6
    for (j in seq_len(ncol(X))) {
      v <- X[, j]
      v[!is.finite(v) | v <= 0] <- x_min[j] / 2
      X[, j] <- v
    }
    cat("Imputation: half minimum per feature.\n")
  } else if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE))
      stop("需要 impute 包用于 KNN 填补。")
    Xt <- t(X)
    im <- impute::impute.knn(data = Xt)$data
    X <- t(im)
    cat("Imputation: KNN.\n")
  }
  X
}

# ---------- 归一化 ----------
#' @keywords internal
ms_normalize <- function(X, method = "median") {
  if (method == "TIC") {
    s <- rowSums(X, na.rm = TRUE)
    s[s == 0] <- 1
    X <- X / s * median(s)
    cat("Normalization: TIC.\n")
  } else if (method == "median") {
    med <- apply(X, 1, median, na.rm = TRUE)
    med[med == 0] <- 1
    X <- X / med * median(med)
    cat("Normalization: sample median.\n")
  } else if (method == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE))
      stop("需要 preprocessCore 包用于分位数归一化。")
    X <- t(preprocessCore::normalize.quantiles(t(X)))
    cat("Normalization: quantile.\n")
  }
  X
}

# ---------- QC 漂移校正 ----------
#' @keywords internal
ms_qc_drift_correction <- function(X, qc_flag, sample_meta) {
  # 假设样本行顺序即注入顺序；若有注入顺序列可替换
  inj_order <- seq_len(nrow(X))
  for (j in seq_len(ncol(X))) {
    y <- X[, j]
    # 仅基于 QC 拟合
    if (sum(qc_flag) >= 5) {
      xx <- inj_order[qc_flag]
      yy <- y[qc_flag]
      # LOESS/样条：用 loess 稳健拟合趋势
      fit <- try(suppressWarnings(stats::loess(yy ~ xx, span = 0.5, degree = 2)), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        trend <- stats::predict(fit, inj_order)
        # 避免负/零
        trend[!is.finite(trend) | trend <= 0] <- median(yy, na.rm = TRUE)
        X[, j] <- y / trend * median(yy, na.rm = TRUE)
      }
    }
  }
  cat("QC-based drift correction: LOESS per feature.\n")
  X
}

# ---------- 基于 QC 的 RSD 过滤 ----------
#' @keywords internal
ms_qc_rsd_filter <- function(X, qc_flag, rsd_cut = 0.30) {
  if (sum(qc_flag) < 3) return(X)
  rsd <- apply(X[qc_flag, , drop = FALSE], 2, function(v) {
    m <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
    if (!is.finite(m) || m == 0) return(1e6)
    s / m
  })
  keep <- rsd <= rsd_cut
  Xf <- X[, keep, drop = FALSE]
  cat(sprintf("QC-RSD filter kept %d/%d features (RSD <= %.2f).\n",
              ncol(Xf), ncol(X), rsd_cut))
  Xf
}

# ========================= 3. 差异分析 =========================

#' @keywords internal
ms_differential <- function(
    X, sample_meta, group_var,
    method = "limma",
    qc_flag_var = NULL,      # 如有列是 TRUE/FALSE 标记 QC，就填列名；否则置 NULL
    qc_label = "QC",         # 如果 QC 在 group_var 中就是一个水平（如 "QC"），写这里
    sort_by = c("none","F")  # limma::topTable 排序方式
) {
  sort_by <- match.arg(sort_by)

  # ---- 0) 识别并剔除 QC 样本 ----
  groups <- as.factor(sample_meta[[group_var]])
  keep_idx <- rep(TRUE, nrow(sample_meta))

  if (!is.null(qc_flag_var) && qc_flag_var %in% colnames(sample_meta)) {
    keep_idx <- keep_idx & !isTRUE(sample_meta[[qc_flag_var]])
  }
  if (!is.null(qc_label) && qc_label %in% levels(groups)) {
    keep_idx <- keep_idx & (groups != qc_label)
  }

  if (!any(keep_idx)) stop("过滤 QC 后没有可用于差异分析的样本。")

  X_use <- X[keep_idx, , drop = FALSE]
  sample_use <- droplevels(sample_meta[keep_idx, , drop = FALSE])
  groups_use <- droplevels(as.factor(sample_use[[group_var]]))

  if (nlevels(groups_use) < 2) {
    stop("过滤 QC 后分组少于 2 个，无法进行差异分析。")
  }

  # ---- 1) limma 或 Wilcoxon ----
  if (method == "limma") {
    if (!requireNamespace("limma", quietly = TRUE))
      stop("需要 limma 包。")

    # 设计矩阵 ~0 + group
    design <- stats::model.matrix(~ 0 + groups_use)
    colnames(design) <- levels(groups_use)

    fit <- limma::lmFit(t(X_use), design)

    if (ncol(design) == 2) {
      # 两组：第二组 - 第一组（按需要可用 relevel 固定参考组）
      contr <- paste0(colnames(design)[2], "-", colnames(design)[1])
      cont <- limma::makeContrasts(contrasts = contr, levels = design)
      fit2 <- limma::contrasts.fit(fit, cont)
      fit2 <- limma::eBayes(fit2)

      tt <- limma::topTable(fit2, number = Inf, sort.by = sort_by)

      out <- data.frame(
        Metabolite     = rownames(tt),
        pvalue         = tt$P.Value,
        padj           = tt$adj.P.Val,
        log2FoldChange = tt$logFC,
        AveExpr        = tt$AveExpr,
        stringsAsFactors = FALSE
      )

    } else {
      # 多组：给出全局模型下各系数结果（不提供 log2FC）
      fit2 <- limma::eBayes(fit)
      tt <- limma::topTable(fit2, number = Inf, sort.by = sort_by)

      out <- data.frame(
        Metabolite = rownames(tt),
        pvalue     = tt$P.Value,
        padj       = tt$adj.P.Val,
        AveExpr    = tt$AveExpr,
        stringsAsFactors = FALSE
      )
    }

  } else { # ---- Wilcoxon：仅两组 ----
    if (nlevels(groups_use) != 2) {
      stop("Wilcoxon 仅支持两组比较。请使用 limma 或指定二分类数据。")
    }
    g1 <- levels(groups_use)[1]; g2 <- levels(groups_use)[2]

    pvs <- apply(X_use, 2, function(v) {
      tryCatch({
        stats::wilcox.test(v[groups_use == g1], v[groups_use == g2])$p.value
      }, error = function(e) NA_real_)
    })
    fc <- apply(X_use, 2, function(v) {
      m1 <- stats::median(v[groups_use == g2], na.rm = TRUE)
      m0 <- stats::median(v[groups_use == g1], na.rm = TRUE)
      log2((m1 + 1e-6) / (m0 + 1e-6))
    })

    out <- data.frame(
      Metabolite     = colnames(X_use),
      pvalue         = pvs,
      padj           = stats::p.adjust(pvs, method = "BH"),
      log2FoldChange = fc,
      AveExpr        = colMeans(X_use, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  out[!is.na(out$pvalue), , drop = FALSE]
}

# ========================= 4. 机器学习重要性 =========================

#' @keywords internal
ms_ml_importance <- function(X, sample_meta, group_var,
                             method = "rf", nfolds = 5, ntree = 500) {
  y <- factor(sample_meta[[group_var]])
  ml_df <- data.frame(X, .Group = y, check.names = FALSE)

  if (method == "rf") {
    if (!requireNamespace("caret", quietly = TRUE))
      stop("需要 caret 包。")
    set.seed(123)
    ctrl <- caret::trainControl(method = "cv", number = nfolds, classProbs = TRUE)
    rf_model <- caret::train(.Group ~ ., data = ml_df, method = "rf",
                             trControl = ctrl, ntree = ntree, importance = TRUE,
                             metric = "Accuracy")
    imp <- caret::varImp(rf_model, scale = TRUE)$importance
    importance_df <- data.frame(
      Metabolite = rownames(imp),
      Importance = imp[,1],
      stringsAsFactors = FALSE
    )
    cat(sprintf("RF CV accuracy: %.2f%%\n", max(rf_model$results$Accuracy) * 100))
  } else if (method == "gbm") {
    if (!requireNamespace("gbm", quietly = TRUE))
      stop("需要 gbm 包。")
    set.seed(123)
    ctrl <- caret::trainControl(method = "cv", number = nfolds, classProbs = TRUE)
    gbm_model <- caret::train(.Group ~ ., data = ml_df, method = "gbm",
                              trControl = ctrl, verbose = FALSE)
    imp <- caret::varImp(gbm_model, scale = TRUE)$importance
    importance_df <- data.frame(
      Metabolite = rownames(imp),
      Importance = imp[,1],
      stringsAsFactors = FALSE
    )
  } else if (method == "boruta") {
    if (!requireNamespace("Boruta", quietly = TRUE))
      stop("需要 Boruta 包。")
    set.seed(123)
    bor <- Boruta::Boruta(.Group ~ ., data = ml_df, doTrace = 0)
    importance_df <- data.frame(
      Metabolite = names(bor$ImpHistory[,1]),
      Importance = apply(bor$ImpHistory, 1, median, na.rm = TRUE),
      Decision  = bor$finalDecision,
      stringsAsFactors = FALSE
    )
  }

  # 标准化 0-1
  if (isTRUE(max(importance_df$Importance) > min(importance_df$Importance))) {
    importance_df$Importance_scaled <- (importance_df$Importance - min(importance_df$Importance)) /
      (max(importance_df$Importance) - min(importance_df$Importance))
  } else {
    importance_df$Importance_scaled <- 0.5
  }
  importance_df
}

# ========================= 5. 相关网络 =========================

#' @keywords internal
ms_network <- function(X, cor_method = "spearman", cor_threshold = 0.6, p_threshold = 0.05) {
  if (requireNamespace("Hmisc", quietly = TRUE)) {
    cm <- Hmisc::rcorr(as.matrix(X), type = cor_method)
    r  <- cm$r; p <- cm$P
  } else {
    r <- cor(X, method = cor_method, use = "pairwise.complete.obs")
    p <- matrix(0.05, nrow = ncol(X), ncol = ncol(X),
                dimnames = list(colnames(X), colnames(X)))
  }
  diag(r) <- 0; diag(p) <- 1
  adj <- (abs(r) >= cor_threshold) & (p < p_threshold)
  adj_num <- adj * 1

  df <- data.frame(
    Metabolite = colnames(X),
    Degree = colSums(adj_num),
    Weighted_degree = colSums(abs(r) * adj_num),
    stringsAsFactors = FALSE
  )

  if (sum(adj_num) > 0 && requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(adj_num * abs(r), mode = "undirected", weighted = TRUE, diag = FALSE)
    df$Betweenness <- igraph::betweenness(g)
    df$Closeness   <- igraph::closeness(g)
    df$Eigenvector <- igraph::eigen_centrality(g)$vector
    df$Hub_score   <- igraph::hub_score(g)$vector
    if (igraph::vcount(g) > 10) {
      comm <- igraph::cluster_fast_greedy(g)
      df$Module <- igraph::membership(comm)
    } else {
      df$Module <- 1
    }
  } else {
    df$Betweenness <- 0; df$Closeness <- 0; df$Eigenvector <- 0; df$Hub_score <- 0; df$Module <- 1
  }

  cat(sprintf("Network edges: %d (|r|>=%.2f, p<%.3f)\n", sum(adj_num)/2, cor_threshold, p_threshold))
  df
}

# ========================= 6. 丰度与流行率 =========================

#' @keywords internal
ms_abundance_metrics <- function(X) {
  # 这里的"丰度"用强度均值，"流行率"用>0 的比例
  data.frame(
    Metabolite = colnames(X),
    Mean_intensity = colMeans(X, na.rm = TRUE),
    Median_intensity = apply(X, 2, median, na.rm = TRUE),
    Max_intensity = apply(X, 2, max, na.rm = TRUE),
    SD_intensity = apply(X, 2, sd, na.rm = TRUE),
    Prevalence = colMeans(X > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# ========================= 7. 综合加权排名 =========================

#' @keywords internal
ms_weighted_ranking <- function(diff_results, importance_df, network_df, abundance_df, weights) {
  # 归一化权重
  s <- sum(unlist(weights)); if (s <= 0) s <- 1
  w <- lapply(weights, function(x) x / s)

  # 合并（Metabolite 作为键）
  z <- merge(diff_results, importance_df, by = "Metabolite", all = TRUE)
  z <- merge(z, network_df, by = "Metabolite", all = TRUE)
  z <- merge(z, abundance_df, by = "Metabolite", all = TRUE)

  num_cols <- sapply(z, is.numeric)
  z[, num_cols][is.na(z[, num_cols])] <- 0

  normalize <- function(x) {
    if (length(unique(x)) <= 1) return(rep(0.5, length(x)))
    rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (rng == 0) return(rep(0.5, length(x)))
    (x - min(x, na.rm = TRUE)) / rng
  }

  # 差异得分
  z$Diff_score_p  <- normalize(1 - z$pvalue)
  z$Diff_score_fc <- normalize(abs(z$log2FoldChange))
  z$Diff_score    <- (z$Diff_score_p + z$Diff_score_fc)/2

  # 重要性
  z$Importance_score <- if ("Importance_scaled" %in% names(z)) z$Importance_scaled else normalize(z$Importance)

  # 网络
  if ("Degree" %in% names(z)) {
    z$Network_score_degree  <- normalize(z$Degree)
    z$Network_score_between <- normalize(z$Betweenness)
    z$Network_score <- (z$Network_score_degree + z$Network_score_between)/2
  } else {
    z$Network_score <- 0.5
  }

  # 丰度（偏好高强度+高出现率）
  if ("Mean_intensity" %in% names(z)) {
    z$Abundance_score_mean <- normalize(z$Mean_intensity)
    z$Abundance_score_prev <- normalize(z$Prevalence)
    z$Abundance_score <- 0.7 * z$Abundance_score_mean + 0.3 * z$Abundance_score_prev
  } else {
    z$Abundance_score <- 0.5
  }

  # 总得分
  z$Weighted_score <-
    w$abundance   * z$Abundance_score +
    w$importance  * z$Importance_score +
    w$differential* z$Diff_score +
    w$network     * z$Network_score

  z <- z[order(z$Weighted_score, decreasing = TRUE), ]
  z$Rank <- seq_len(nrow(z))

  key <- c("Rank","Metabolite","Weighted_score",
           "Mean_intensity","Prevalence",
           "pvalue","padj","log2FoldChange",
           "Importance","Degree","Betweenness")
  key <- key[key %in% names(z)]
  summary_tbl <- z[, key, drop = FALSE]

  list(full_results = z, summary_results = summary_tbl)
}

# ========================= 8. 可视化（可选） =========================

#' Visualize MS Feature Selection Results
#' @param results ms_feature_selection 的返回
#' @param top_n 展示前 N
#' @export
plot_ms_results <- function(results, top_n = 20) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("需要 ggplot2 包。")

  df <- results$full_results
  top_df <- head(df, top_n)

  p1 <- ggplot2::ggplot(top_df,
                        ggplot2::aes(x = stats::reorder(Metabolite, Weighted_score),
                                     y = Weighted_score)) +
    ggplot2::geom_bar(stat = "identity", alpha = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Top", top_n, "Metabolite Importance (Integrated)"),
                  x = "Metabolite", y = "Weighted score")

  p1
}

# ========================= 9. 包信息和方法 =========================

#' Print Method for ms_feature_selection Objects
#'
#' @param x Object of class ms_feature_selection
#' @param ... Additional arguments (not used)
#'
#' @export
print.ms_feature_selection <- function(x, ...) {
  cat("MS Feature Selection Results\n")
  cat("============================\n\n")

  # Analysis parameters
  params <- x$parameters
  cat("Analysis Parameters:\n")
  cat(sprintf("  Group variable: %s\n", params$group_var))
  cat(sprintf("  Differential method: %s\n", params$diff_method))
  cat(sprintf("  ML method: %s\n", params$ml_method))
  cat(sprintf("  Normalization: %s\n", params$normalize_method))
  cat(sprintf("  Correlation method: %s\n", params$cor_method))

  # Data summary
  cat("\nData Summary:\n")
  cat(sprintf("  Total features analyzed: %d\n", nrow(x$full_results)))
  cat(sprintf("  Features in final ranking: %d\n", nrow(x$final_ranking)))

  # Results summary
  if ("pvalue" %in% names(x$diff_results)) {
    sig_features <- sum(x$diff_results$pvalue < 0.05, na.rm = TRUE)
    cat(sprintf("  Significantly different features (p < 0.05): %d\n", sig_features))
  }

  if ("Degree" %in% names(x$network_results)) {
    connected_features <- sum(x$network_results$Degree > 0)
    cat(sprintf("  Network-connected features: %d\n", connected_features))
  }

  # Top features
  cat("\nTop 5 Features (by weighted score):\n")
  top5 <- head(x$final_ranking, 5)
  for (i in 1:min(5, nrow(top5))) {
    cat(sprintf("  %d. %s (score: %.3f)\n",
                i, top5$Metabolite[i], top5$Weighted_score[i]))
  }

  cat("\nUse plot_ms_results() to visualize results\n")
}

#' Summary Method for ms_feature_selection Objects
#'
#' @param object Object of class ms_feature_selection
#' @param ... Additional arguments (not used)
#'
#' @export
summary.ms_feature_selection <- function(object, ...) {
  print(object)

  # Additional detailed statistics
  cat("\nDetailed Statistics:\n")
  cat("===================\n")

  # Component score distributions
  full_res <- object$full_results
  if ("Diff_score" %in% names(full_res)) {
    cat(sprintf("Differential scores - Mean: %.3f, SD: %.3f\n",
                mean(full_res$Diff_score, na.rm = TRUE),
                sd(full_res$Diff_score, na.rm = TRUE)))
  }

  if ("Importance_score" %in% names(full_res)) {
    cat(sprintf("Importance scores - Mean: %.3f, SD: %.3f\n",
                mean(full_res$Importance_score, na.rm = TRUE),
                sd(full_res$Importance_score, na.rm = TRUE)))
  }

  if ("Network_score" %in% names(full_res)) {
    cat(sprintf("Network scores - Mean: %.3f, SD: %.3f\n",
                mean(full_res$Network_score, na.rm = TRUE),
                sd(full_res$Network_score, na.rm = TRUE)))
  }

  if ("Abundance_score" %in% names(full_res)) {
    cat(sprintf("Abundance scores - Mean: %.3f, SD: %.3f\n",
                mean(full_res$Abundance_score, na.rm = TRUE),
                sd(full_res$Abundance_score, na.rm = TRUE)))
  }
}
