


.normalize_imp <- function(imp, method = c("l1","softmax","minmax","none")) {
  method <- match.arg(method)
  if (method == "none") return(imp)

  if (method == "l1") {
    imp_pos <- pmax(imp, 0)                     # 先把负数截为 0，避免负贡献
    s <- sum(imp_pos, na.rm = TRUE)
    if (s > 0) {
      return(imp_pos / s)                       # 和为 1
    } else {
      # 如果全是 0 或 NA，则退化为 softmax，避免全 0
      w <- exp(imp - max(imp, na.rm = TRUE))
      return(w / sum(w, na.rm = TRUE))
    }
  }

  if (method == "softmax") {
    w <- exp(imp - max(imp, na.rm = TRUE))
    return(w / sum(w, na.rm = TRUE))
  }

  if (method == "minmax") {
    mn <- min(imp, na.rm = TRUE); mx <- max(imp, na.rm = TRUE)
    if (mx > mn) {
      return((imp - mn) / (mx - mn))
    } else {
      return(rep(0, length(imp)))
    }
  }
}

#' @title Random Forest-based Microbial Association Network
#'
#' @description
#' Construct a microbial association network using Random Forests from a
#' `phyloseq` object. Each taxon (OTU/ASV/species) is iteratively treated as the
#' response variable, with the remaining taxa as predictors. Variable
#' importance is extracted from the Random Forest model, normalized such that
#' the importance scores for all predictors sum to 1, and optionally evaluated
#' for significance using permutation testing.
#'
#' @param ps A `phyloseq` object containing OTU/ASV table.
#' @param ntree Integer. Number of trees to grow in each Random Forest.
#'   Default is \code{500}.
#' @param normalize Logical. Whether to normalize variable importance so that
#'   scores sum to 1 for each model. Default is \code{TRUE}.
#' @param scale.method Character. Method used for normalization.
#'   Options: \code{"l1"}, \code{"softmax"}, \code{"minmax"}, \code{"none"}.
#'   Default is \code{"l1"}.
#' @param n_perm Integer. Number of permutations for significance testing.
#'   If > 0, permutation-based p-values are computed.
#'   If 0, no permutation test is performed. Default is \code{100}.
#' @param p.threshold Numeric. Significance threshold for permutation test.
#'   Edges with p-value < threshold are retained. Default is \code{0.05}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{123}.
#'
#' @details
#' For each taxon, a Random Forest model is trained using all other taxa as
#' predictors. Variable importance scores are extracted using permutation
#' importance. Importance scores can be normalized by several methods:
#' \itemize{
#'   \item \code{"l1"}: Divide by L1 norm (absolute sum), making scores sum to 1.
#'   \item \code{"softmax"}: Convert to probabilities via softmax.
#'   \item \code{"minmax"}: Rescale to [0, 1] using min-max normalization.
#'   \item \code{"none"}: Use raw importance values.
#' }
#'
#' When \code{n_perm > 0}, a permutation test is performed by randomly shuffling
#' the response variable \code{y} and recomputing importance values across
#' \code{n_perm} iterations. P-values are then computed as the proportion of
#' null importance values greater than or equal to the observed importance.
#'
#' @return A list with two components:
#' \describe{
#'   \item{edges}{A data frame of significant associations with columns:
#'     \code{target}, \code{predictor}, \code{importance}, and (if applicable) \code{pval}.}
#'   \item{importance_all}{A data frame of all importance values (before significance filtering).}
#' }
#'
#' @author Wen Tao <\email{2018203048@@njau.edu.cn}>
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("GlobalPatterns")
#'
#' # Subset to 20 taxa for demonstration
#' ps_sub <- prune_taxa(taxa_names(GlobalPatterns)[1:20], GlobalPatterns)
#'
#' res <- rf_network_from_phyloseq(
#'   ps = ps_sub,
#'   ntree = 100,
#'   normalize = TRUE,
#'   scale.method = "l1",
#'   n_perm = 10,
#'   p.threshold = 0.05
#' )
#'
#' head(res$edges)
#' head(res$importance_all)
#' }
#'
#' @export
rf_network_from_phyloseq <- function(ps,
                                     ntree = 500,
                                     normalize = TRUE,          # 是否归一化
                                     scale.method = c("l1","softmax","minmax","none"),
                                     n_perm = 100,             # 置换次数；>0 则做显著性检验
                                     p.threshold = 0.05,
                                     seed = 123) {
  set.seed(seed)
  scale.method <- match.arg(scale.method)

  # 1) 取 OTU 表（样本在行，特征在列）
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- as.data.frame(otu, check.names = FALSE)

  taxa_names_vec <- colnames(otu)
  n_taxa <- ncol(otu)

  edge_list <- vector("list", n_taxa)
  imp_list  <- vector("list", n_taxa)

  # 2) 循环每个微生物作为 y
  for (i in seq_len(n_taxa)) {
    y_name <- taxa_names_vec[i]
    y <- otu[[y_name]]
    X <- otu[, -i, drop = FALSE]
    all_pred <- colnames(X)

    df <- data.frame(y = y, X, check.names = FALSE)

    # 3) 训练 RF（用列名 "y" 作为因变量，避免特殊字符问题）
    rf <- ranger(
      dependent.variable.name = "y",
      data = df,
      num.trees = ntree,
      importance = "permutation"
    )
    imp_raw <- rf$variable.importance
    # 对齐所有预测变量（若个别变量被模型丢弃，补 0）
    imp_raw <- imp_raw[all_pred]; imp_raw[is.na(imp_raw)] <- 0

    # 4) 归一化（使当前 y 下所有 importance 之和为 1）
    imp_use <- if (isTRUE(normalize)) .normalize_imp(imp_raw, scale.method) else imp_raw

    imp_df <- data.frame(
      target = y_name,
      predictor = all_pred,
      importance = as.numeric(imp_use),
      stringsAsFactors = FALSE
    )

    # 5) 显著性：置换检验（对每次置换也做同样的归一化）
    if (!is.null(n_perm) && n_perm > 0) {
      perm_mat <- replicate(n_perm, {
        df$y <- sample(df$y)  # 打乱 y
        rf_p <- ranger(
          dependent.variable.name = "y",
          data = df,
          num.trees = ntree,
          importance = "permutation"
        )
        imp_p <- rf_p$variable.importance
        imp_p <- imp_p[all_pred]; imp_p[is.na(imp_p)] <- 0
        if (isTRUE(normalize)) imp_p <- .normalize_imp(imp_p, scale.method)
        imp_p
      })
      # 确保是矩阵（行=变量，列=置换）
      if (is.null(dim(perm_mat))) perm_mat <- matrix(perm_mat, nrow = length(all_pred))
      rownames(perm_mat) <- all_pred

      # p 值：null >= real 的比例
      pvals <- vapply(all_pred, function(v) {
        real <- imp_use[v]
        nulls <- perm_mat[v, ]
        mean(nulls >= real, na.rm = TRUE)
      }, numeric(1))

      imp_df$pval <- pvals
      sig_df <- dplyr::filter(imp_df, pval < p.threshold)

    } else {
      # 不做置换：按分位数筛选（在归一化后的 importance 上）
      thr <- stats::quantile(imp_use, probs = 1 - select_top_p, na.rm = TRUE)
      sig_df <- dplyr::filter(imp_df, importance >= thr)
    }

    edge_list[[i]] <- sig_df
    imp_list[[i]]  <- imp_df
  }

  edges <- dplyr::bind_rows(edge_list)
  importance_all <- dplyr::bind_rows(imp_list)

  # 输出：边表 + 全量重要性
  list(
    edges = edges,                 # 列：target, predictor, importance(和为1), 可含 pval
    importance_all = importance_all
  )
}
