#' MicrobeFinderFromMetabolite: find microbes associated with a target metabolite
#'
#' 给定一个代谢物（非靶向代谢组矩阵中的一列），在 phyloseq 微生物数据中寻找潜在互作的微生物特征。
#' 采用多维证据：相关（CLR/SparCC 可选）、网络中心性、差异丰度（代谢物高低分组）、机器学习
#' （用微生物预测该代谢物，得到微生物重要性），以及丰度/流行率优先。
#'
#' 依赖包可选：DESeq2（差异分析更稳），mixOmics（可做 PLS 回归 VIP），
#' igraph（网络）、glmnet / randomForest（机器学习）。没装会自动降级。
#'
#' @param phyloseq_obj phyloseq对象（微生物计数或相对丰度；建议原始计数）
#' @param metabolome_mat 矩阵/数据框（样本×代谢物）
#' @param target_metabolite_name 字符串：目标代谢物列名
#' @param top_n 返回Top微生物数
#' @param cor_method 相关方法："spearman"（默认，基于CLR）、"pearson"、"sparcc"（需 SpiecEasi）
#' @param prevalence_threshold 微生物流行率阈值（默认0.1）
#' @param detection_threshold 相对丰度检出阈值（默认1e-4；用于计算流行率）
#' @param use_clr 是否对微生物做CLR（默认TRUE；DESeq2差异仍用原始计数）
#' @param group_split 策略将代谢物分为High/Low：c("median","tertile","quantile_0.4_0.6")
#' @return list：top_microbes, all_scores, correlation_results, network_metrics,
#'               diff_results, ml_results, aligned_data（样本对齐/变换后的数据）
#' @export
#'
    # group_split ="median"
find_associated_microbes <- function(
    phyloseq_obj,
    metabolome_mat,
    target_metabolite_name,
    top_n = 30,
    cor_method = c("spearman","pearson","sparcc"),
    prevalence_threshold = 0.10,
    detection_threshold  = 1e-4,
    use_clr = TRUE,
    group_split = c("median","tertile","quantile_0.4_0.6")
){
  cor_method  <- match.arg(cor_method)
  group_split <- match.arg(group_split)

  # ---------- 1) 对齐样本 & 基础准备 ----------
  dat <- .prepare_data_micro_side(phyloseq_obj, metabolome_mat, target_metabolite_name,
                                  prevalence_threshold, detection_threshold, use_clr)
  X_clr  <- dat$micro_clr        # 样本×微生物（CLR后；做相关和ML）
  X_rel  <- dat$micro_rel        # 样本×微生物（相对丰度；做流行率/网络）
  X_cnt  <- dat$micro_counts     # 微生物原始计数（做DESeq2）
  y_vec  <- dat$metab_vec        # 目标代谢物（连续）
  samples <- dat$samples
  message(sprintf("Aligned samples: %d; microbes: %d; metabolite: %s",
                  length(samples), ncol(X_clr), target_metabolite_name))

  # ---------- 2) 相关分析 ----------
  cor_res <- .cor_micro_vs_metab(X = X_clr, y = y_vec, method = cor_method)

  # ---------- 3) 微生物共现网络（基于相对丰度或CLR） ----------
  net_res <- .network_microbes(X_rel)  # 用相对丰度更稳；也可改 X_clr

  # ---------- 4) 差异分析（以代谢物高/低组作为分组） ----------
  grp <- .split_by_metabolite(y_vec, strategy = group_split)
  diff_res <- .diff_micro_by_group(phyloseq_obj, samples, grp)  # 优先 DESeq2; fallback Wilcoxon

  # ---------- 5) 机器学习：用微生物（X_clr）预测代谢物（y） ----------
  ml_res <- .ml_importance_micro_to_metab(X = X_clr, y = y_vec)

  # ---------- 6) 丰度/流行率指标 ----------
  abn_res <- .abundance_prevalence(X_rel)

  # ---------- 7) 多指标集成打分（与您现有 integrated_scoring 风格一致） ----------
  scores <- .integrated_scoring_micro(
    cor_results = cor_res,
    network_metrics = net_res,
    diff_results = diff_res,
    ml_results = ml_res,
    abn_results = abn_res,
    X_rel = X_rel
  )

  top_tbl <- head(scores, top_n)

  out <- list(
    top_microbes        = top_tbl,
    all_scores          = scores,
    correlation_results = cor_res,
    network_metrics     = net_res,
    diff_results        = diff_res,
    ml_results          = ml_res,
    aligned_data        = list(X_clr = X_clr, X_rel = X_rel, X_counts = X_cnt,
                               y = y_vec, samples = samples)
  )
  class(out) <- "microbe_association_results"
  return(out)
}

# -------------------- 内部函数 --------------------

# 对齐 & 预处理（过滤低流行率、做CLR、保留相对丰度与计数）
.prepare_data_micro_side <- function(phyloseq_obj, metabolome_mat, metab_name = target_metabolite_name,
                                     prevalence_threshold, detection_threshold, use_clr){
  otu <- as.data.frame(phyloseq::otu_table(phyloseq_obj))
  taxa_are_rows <- phyloseq::taxa_are_rows(phyloseq_obj)
  if (taxa_are_rows) otu <- t(otu) # 样本×微生物

  # 代谢物矩阵：样本×代谢物
  M <- if (is.data.frame(metabolome_mat)) as.matrix(metabolome_mat) else metabolome_mat
  if (!is.matrix(M)) M <- as.matrix(M)

  common_samples <- intersect(rownames(otu), rownames(M))
  if (length(common_samples) == 0) stop("No common samples between phyloseq and metabolome_mat.")

  otu <- otu[common_samples, , drop = FALSE]
  M   <- M[common_samples, , drop = FALSE]

  if (!metab_name %in% colnames(M)) stop("Target metabolite not found in metabolome_mat.")
  y <- as.numeric(M[, metab_name])

  # 相对丰度（TSS）
  otu_rel <- sweep(otu, 1, pmax(rowSums(otu), 1), "/")

  # 流行率过滤
  prev <- colMeans(otu_rel > detection_threshold, na.rm = TRUE)
  keep <- prev >= prevalence_threshold
  if (sum(keep) < 5) {
    warning("After prevalence filtering, <5 taxa remain; relaxing threshold.")
    keep <- order(prev, decreasing = TRUE)[1:min(50, ncol(otu_rel))]
  }
  otu_rel <- otu_rel[, keep, drop = FALSE]
  otu_cnt <- as.data.frame(otu)[, keep, drop = FALSE]

  # CLR（用于相关&ML）
  if (use_clr) {
    psd <- apply(otu_rel, 2, function(z) any(z > 0))
    otu_rel <- otu_rel[, psd, drop = FALSE]
    clr <- t(apply(otu_rel, 1, function(v){
      v <- v + 1e-6
      log(v) - mean(log(v))
    }))
  } else {
    clr <- otu_rel
  }

  # 同步计数矩阵列
  otu_cnt <- otu_cnt[, colnames(otu_rel), drop = FALSE]

  list(
    micro_rel   = otu_rel,
    micro_clr   = clr,
    micro_counts= otu_cnt,
    metab_vec   = y,
    samples     = common_samples
  )
}

# 相关：Spearman/Pearson（SparCC 可选）
.cor_micro_vs_metab <- function(X, y, method = "spearman"){
  p <- ncol(X)
  res <- data.frame(
    microbe = colnames(X),
    cor     = NA_real_,
    pval    = NA_real_,
    stringsAsFactors = FALSE
  )

  if (method %in% c("spearman","pearson")){
    for (j in seq_len(p)){
      ct <- suppressWarnings(cor.test(X[, j], y, method = method))
      res$cor[j]  <- unname(ct$estimate)
      res$pval[j] <- ct$p.value
    }
  } else if (method == "sparcc"){
    if (!requireNamespace("SpiecEasi", quietly = TRUE)){
      warning("SpiecEasi not installed; fallback to Spearman.")
      return(.cor_micro_vs_metab(X, y, method = "spearman"))
    }
    # 对每个特征与 y 的相关，SparCC 不是单变量对连续 y 的；这里折中：仍用 spearman
    for (j in seq_len(p)){
      ct <- suppressWarnings(cor.test(X[, j], y, method = "spearman"))
      res$cor[j]  <- unname(ct$estimate)
      res$pval[j] <- ct$p.value
    }
  }
  res$fdr <- p.adjust(res$pval, method = "BH")
  res
}

# 微生物网络（共现）：用相关阈值构网并算中心性
.network_microbes <- function(X_rel, cor_threshold = 0.3){
  C <- suppressWarnings(cor(X_rel, use = "pairwise.complete.obs", method = "spearman"))
  C[!is.finite(C)] <- 0
  diag(C) <- 0
  C[abs(C) < cor_threshold] <- 0

  if (!requireNamespace("igraph", quietly = TRUE)){
    warning("igraph not installed; network metrics set to zero.")
    return(data.frame(
      microbe     = colnames(X_rel),
      degree      = 0, betweenness = 0, closeness = 0,
      eigenvector = 0, hub_score  = 0,
      stringsAsFactors = FALSE
    ))
  }

  g <- igraph::graph_from_adjacency_matrix(abs(C), mode = "undirected", weighted = TRUE)
  data.frame(
    microbe     = colnames(X_rel),
    degree      = igraph::degree(g),
    betweenness = igraph::betweenness(g, weights = NA),
    closeness   = igraph::closeness(g, weights = NA),
    eigenvector = igraph::eigen_centrality(g, weights = NA)$vector,
    hub_score   = igraph::hub_score(g, weights = NA)$vector,
    stringsAsFactors = FALSE
  )
}

# 以代谢物高/低组为设计做差异丰度：优先 DESeq2（原始计数），没有就 Wilcoxon（相对丰度）
# phyloseq_obj
# keep_samples =  samples
# groups = grp


.diff_micro_by_group <- function(phyloseq_obj, keep_samples, groups){
  # groups: factor, length = |keep_samples|
  ps <- phyloseq::prune_samples(keep_samples, phyloseq_obj)
  md <- as.data.frame(phyloseq::sample_data(ps))
  md$.metab_group <- factor(groups)
  phyloseq::sample_data(ps) <- phyloseq::sample_data(md)

  if (requireNamespace("DESeq2", quietly = TRUE)){
    dds <- phyloseq::phyloseq_to_deseq2(ps, ~ .metab_group)
    # sizeFactors
    geoMeans <- apply(DESeq2::counts(dds), 1, function(row) {
      if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
    })
    dds <- DESeq2::estimateSizeFactors(dds, geoMeans = geoMeans)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    resultsNames(dds)
    res <- DESeq2::results(dds, name = ".metab_group_Low_vs_High")
    out <- data.frame(
      microbe = rownames(res),
      pvalue  = res$pvalue,
      padj    = res$padj,
      log2FC  = res$log2FoldChange,
      baseMean= res$baseMean,
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$pvalue), ]
    return(out)
  }

  # fallback: Wilcoxon on relative abundance
  otu <- as.data.frame(phyloseq::otu_table(ps))
  if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  otu_rel <- sweep(otu, 1, pmax(rowSums(otu), 1), "/")
  grp <- factor(groups)
  g1 <- levels(grp)[1]; g2 <- levels(grp)[2]
  pvs <- apply(otu_rel, 2, function(v){
    tryCatch(stats::wilcox.test(v[grp==g1], v[grp==g2])$p.value, error = function(e) NA_real_)
  })
  fc <- apply(otu_rel, 2, function(v){
    m1 <- median(v[grp==g2], na.rm = TRUE); m0 <- median(v[grp==g1], na.rm = TRUE)
    log2((m1 + 1e-6)/(m0 + 1e-6))
  })
  data.frame(
    microbe = colnames(otu_rel),
    pvalue  = pvs,
    padj    = p.adjust(pvs, method = "BH"),
    log2FC  = fc,
    baseMean= colMeans(otu_rel),
    stringsAsFactors = FALSE
  )
}

# 将连续代谢物分组 High/Low
.split_by_metabolite <- function(y, strategy = c("median","tertile","quantile_0.4_0.6")){
  strategy <- match.arg(strategy)
  if (strategy == "median"){
    thr <- median(y, na.rm = TRUE)
    factor(ifelse(y > thr, "High","Low"))
  } else if (strategy == "tertile"){
    q <- quantile(y, probs = c(1/3, 2/3), na.rm = TRUE)
    cut(y, breaks = c(-Inf, q[1], q[2], Inf), labels = c("Low","Mid","High")) %>%
      droplevels() %>% forcats::fct_collapse(High = c("High"), Low = c("Low"))
  } else {
    q <- quantile(y, probs = c(0.4, 0.6), na.rm = TRUE)
    zz <- y
    zz[y <= q[1] | y >= q[2]] -> sel
    # 仅保留两端；中间的标 NA
    grp <- rep(NA_character_, length(y))
    grp[y <= q[1]] <- "Low"
    grp[y >= q[2]] <- "High"
    factor(grp)
  }
}

# 机器学习：用微生物预测代谢物，输出微生物重要性
.ml_importance_micro_to_metab <- function(X, y){
  out <- list()

  # 1) RandomForest 回归
  if (requireNamespace("randomForest", quietly = TRUE)){
    set.seed(123)
    rf <- randomForest::randomForest(x = X, y = y, ntree = 500, importance = TRUE)
    imp <- as.data.frame(randomForest::importance(rf))
    imp$microbe <- rownames(imp)
    out$rf_importance <- imp[, c("microbe","IncNodePurity"), drop = FALSE]
    colnames(out$rf_importance)[2] <- "rf_importance"
  } else {
    out$rf_importance <- NULL
  }

  # 2) LASSO 回归
  if (requireNamespace("glmnet", quietly = TRUE)){
    cv <- glmnet::cv.glmnet(as.matrix(X), y, alpha = 1, family = "gaussian")
    coefm <- coef(cv, s = "lambda.min")
    lasso <- data.frame(microbe = rownames(coefm)[-1],
                        coef = as.numeric(coefm[-1, ]))
    out$lasso_coef <- lasso
  } else {
    out$lasso_coef <- NULL
  }

  # 3) PLS 回归 VIP（可选）
  if (requireNamespace("mixOmics", quietly = TRUE)){
    pls <- tryCatch(mixOmics::pls(X, y, ncomp = 2), error = function(e) NULL)
    if (!is.null(pls)) {
      vip <- tryCatch(mixOmics::vip(pls), error = function(e) NULL)
      if (!is.null(vip)) {
        out$pls_vip <- data.frame(microbe = rownames(vip),
                                  vip1 = vip[,1], vip2 = if (ncol(vip)>=2) vip[,2] else NA_real_)
      }
    }
  }

  out
}

# 丰度与流行率
.abundance_prevalence <- function(X_rel){
  data.frame(
    microbe    = colnames(X_rel),
    mean_abn   = colMeans(X_rel, na.rm = TRUE),
    prevalence = colMeans(X_rel > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# 集成打分（镜像你上一版 integrated_scoring 的思路）
.integrated_scoring_micro <- function(cor_results, network_metrics, diff_results,
                                      ml_results, abn_results, X_rel,
                                      weights = c(
                                        abundance    = 0.15,
                                        correlation  = 0.25,
                                        network      = 0.15,
                                        differential = 0.20,
                                        ml           = 0.25
                                      )){
  sc <- data.frame(microbe = colnames(X_rel), stringsAsFactors = FALSE)

  # 1) Abundance/Prevalence
  sc <- merge(sc, abn_results, by = "microbe", all.x = TRUE)
  sc$abundance_score <- rank(sc$mean_abn, na.last = "keep")/nrow(sc) * 0.7 +
    rank(sc$prevalence, na.last = "keep")/nrow(sc) * 0.3
  sc$abundance_score[is.na(sc$abundance_score)] <- 0

  # 2) Correlation（绝对值越大越好）
  sc <- merge(sc, cor_results[, c("microbe","cor","fdr")], by = "microbe", all.x = TRUE)
  sc$correlation_score <- abs(sc$cor); sc$correlation_score[is.na(sc$correlation_score)] <- 0

  # 3) Network
  sc <- merge(sc, network_metrics[, c("microbe","degree","betweenness")], by = "microbe", all.x = TRUE)
  sc$network_score <- (rank(sc$degree, na.last = "keep") +
                         rank(sc$betweenness, na.last = "keep"))/(2*nrow(sc))
  sc$network_score[is.na(sc$network_score)] <- 0

  # 4) Differential
  if (!is.null(diff_results) && nrow(diff_results)){
    tmp <- diff_results[, c("microbe","padj","log2FC")]
    sc <- merge(sc, tmp, by = "microbe", all.x = TRUE)
    sc$diff_raw <- -log10(sc$padj + 1e-12) * abs(sc$log2FC)
    sc$diff_score <- rank(sc$diff_raw, na.last = "keep")/nrow(sc)
    sc$diff_score[is.na(sc$diff_score)] <- 0
  } else {
    sc$diff_score <- 0
  }

  # 5) ML
  sc$ml_score <- 0
  if (!is.null(ml_results$rf_importance)){
    sc <- merge(sc, ml_results$rf_importance, by = "microbe", all.x = TRUE)
    sc$ml_score <- rank(sc$rf_importance, na.last = "keep")/nrow(sc)
    sc$ml_score[is.na(sc$ml_score)] <- 0
  }
  if (!is.null(ml_results$lasso_coef)){
    sc <- merge(sc, ml_results$lasso_coef, by = "microbe", all.x = TRUE)
    co <- abs(sc$coef); co[is.na(co)] <- 0
    sc$ml_score <- pmax(sc$ml_score, rank(co)/nrow(sc))
  }
  if (!is.null(ml_results$pls_vip)){
    sc <- merge(sc, ml_results$pls_vip[,c("microbe","vip1")], by = "microbe", all.x = TRUE)
    v <- sc$vip1; v[is.na(v)] <- 0
    sc$ml_score <- pmax(sc$ml_score, rank(v)/nrow(sc))
  }

  # 总分
  sc$integrated_score <-
    weights["abundance"]    * sc$abundance_score +
    weights["correlation"]  * sc$correlation_score +
    weights["network"]      * sc$network_score +
    weights["differential"] * sc$diff_score +
    weights["ml"]           * sc$ml_score

  sc <- sc[order(sc$integrated_score, decreasing = TRUE), ]
  rownames(sc) <- NULL
  sc
}
