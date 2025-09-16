#' Songbird-like multinomial regression in pure R
#'
#' @description
#' This function mimics Songbird (QIIME2, Python) using glmnet multinomial regression.
#' It outputs feature coefficients (log-fold change estimates). Optionally, a "pseudo_p"
#' score is added based on coefficient ranking.
#'
#' @param ps A phyloseq object
#' @param group Column name in sample_data for grouping variable
#' @param alpha Elastic-net mixing parameter (default 0 = ridge, 1 = lasso)
#' @param lambda Regularization parameter (default = lambda.min chosen by CV)
#' @param pseudo_p Logical; if TRUE, output a pseudo p-value based on rank (default FALSE)
#'
#' @return A list with:
#' \itemize{
#'   \item tab.songbirdR: tidy coefficient table (all features, all groups)
#'   \item diff.tab: standardized output (OTU, Group, method, coef, adjust.p/pseudo_p)
#' }
#' @export
#'
#' @examples
#' \dontrun{
#'   library(phyloseq)
#'   data(GlobalPatterns)
#'   # subset a small dataset
#'   ps.sub <- prune_samples(sample_data(GlobalPatterns)$SampleType %in% c("Feces","Skin"), GlobalPatterns)
#'   ps.sub <- prune_taxa(taxa_sums(ps.sub) > 100, ps.sub)
#'
#'   res <- songbirdR.micro(ps.sub, group="SampleType", alpha=0, pseudo_p=TRUE)
#'   head(res$diff.tab)
#' }
songbirdR.micro <- function(ps, group = "Group", alpha = 0, lambda = NULL,
                            pseudo_p = TRUE, collapse = TRUE, sig_alpha = 0.05) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Need glmnet")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Need tidyr")

  otu <- phyloseq::otu_table(ps)
  if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  otu <- as(otu, "matrix")
  meta <- as(phyloseq::sample_data(ps), "data.frame")

  common.samples <- intersect(rownames(otu), rownames(meta))
  otu <- otu[common.samples, , drop = FALSE]
  meta <- meta[common.samples, , drop = FALSE]

  y <- factor(meta[[group]])
  if (nlevels(y) < 2) stop("Grouping variable must have >= 2 levels")

  otu.rel <- sweep(otu, 1, rowSums(otu), "/")
  X <- log(otu.rel + 1e-6)

  fit <- glmnet::cv.glmnet(X, y, family = "multinomial", alpha = alpha)
  if (is.null(lambda)) lambda <- fit$lambda.min
  coef_list <- glmnet::coef.glmnet(fit$glmnet.fit, s = lambda)

  coef_long <- lapply(names(coef_list), function(cls) {
    cf <- as.matrix(coef_list[[cls]])
    data.frame(OTU = rownames(cf),
               Group = cls,
               coef = as.numeric(cf),
               stringsAsFactors = FALSE)
  }) %>% dplyr::bind_rows()

  coef_long <- coef_long %>% dplyr::filter(OTU != "(Intercept)")

  # collapse: 每个 OTU 只保留最大绝对值
  if (collapse) {
    coef_long <- coef_long %>%
      dplyr::group_by(OTU) %>%
      dplyr::slice_max(order_by = abs(coef), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  }

  if (pseudo_p) {
    coef_long <- coef_long %>%
      dplyr::mutate(pseudo_p = rank(-abs(coef)) / n())
    diff.tab <- coef_long %>%
      dplyr::arrange(dplyr::desc(abs(coef))) %>%
      dplyr::mutate(method = "songbirdR",
                    significant = pseudo_p <= sig_alpha) %>%
      dplyr::select(OTU, Group, method, coef, adjust.p = pseudo_p, significant)
  } else {
    diff.tab <- coef_long %>%
      dplyr::arrange(dplyr::desc(abs(coef))) %>%
      dplyr::mutate(method = "songbirdR",
                    adjust.p = NA,
                    significant = FALSE) %>%
      dplyr::select(OTU, Group, method, coef, adjust.p, significant)
  }

  # 只返回显著的
  diff.sig <- diff.tab %>% dplyr::filter(significant)

  return(list(tab.songbirdR = coef_long,
              diff.tab = diff.tab,
              diff.sig = diff.sig))
}

