#' PCA loading analysis for major microbial features
#'
#' Performs PCA on transformed microbiome abundance data, calculates weighted
#' feature importance from PCA loadings, summarizes the abundance of major
#' features, and returns ordination and feature-level plots and tables.
#'
#' @param ps A phyloseq object.
#' @param Top Integer. Number of top PCA-important features retained for display.
#' @param group Character. Sample metadata column containing group labels.
#' @param taxrank Character or NULL. Taxonomic rank used to label features.
#' @param cumulative Numeric. Cumulative explained-variance threshold used to
#'   determine how many PCs contribute to weighted feature importance.
#' @param transform Character. Transformation used before PCA. `"hellinger"`
#'   applies a Hellinger transformation; other values use relative abundance.
#'
#' @return A list containing a combined feature plot, the full feature table,
#'   PCA ordination, scree and importance plots, abundance plots and tables,
#'   group-test results, sample PCA scores, explained variance, loadings,
#'   the `prcomp` object, and the number of PCs used for importance.
#'
#' @family biomarker analysis
#' @export
loadingPCA.metm2 <- function(
    ps, Top = 20, group = "Group", taxrank = NULL,
    cumulative = 0.80, transform = "hellinger"
) {
  ## ---------- data ----------
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)

  meta <- data.frame(
    phyloseq::sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!group %in% names(meta))
    stop("sample_data中不存在: ", group)

  grp <- setNames(
    as.character(meta[[group]]),
    rownames(meta)
  )
  ## 原始相对丰度(%)
  rel <- sweep(otu, 1, rowSums(otu), "/") * 100
  rel[!is.finite(rel)] <- 0

  ## PCA推荐使用Hellinger；展示仍用原始相对丰度
  x <- if (transform == "hellinger") sqrt(rel / 100) else rel / 100
  x <- x[, apply(x, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

  ## ---------- PCA ----------
  pca <- stats::prcomp(x, center = TRUE, scale. = FALSE)
  eig <- pca$sdev^2
  var_exp <- eig / sum(eig)
  cum_exp <- cumsum(var_exp)

  k <- which(cum_exp >= cumulative)[1]
  if (is.na(k)) k <- length(var_exp)
  k <- min(max(k, 2), ncol(pca$rotation))

  ## ---------- feature importance ----------
  load <- pca$rotation
  imp <- rowSums(sweep(load[, 1:k, drop = FALSE]^2, 2, var_exp[1:k], "*"))
  imp <- imp / sum(imp) * 100

  feature <- data.frame(
    TaxonID = rownames(load),
    PC1_loading = load[, 1],
    PC2_loading = load[, 2],
    PC1_contribution = load[, 1]^2 / sum(load[, 1]^2) * 100,
    PC2_contribution = load[, 2]^2 / sum(load[, 2]^2) * 100,
    PCA_importance = imp,
    Mean_abundance = colMeans(rel)[rownames(load)],
    Prevalence = colMeans(otu > 0)[rownames(load)] * 100,
    stringsAsFactors = FALSE
  )

  ## taxonomy名称
  tax <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
  feature$Taxon <- feature$TaxonID

  if (!is.null(tax) && !is.null(taxrank) &&
      taxrank %in% colnames(tax)) {
    tx <- as.data.frame(tax)
    nm <- as.character(tx[feature$TaxonID, taxrank])
    ok <- !is.na(nm) & nzchar(nm) &
      !tolower(nm) %in% c("unassigned", "unknown", "na")
    feature$Taxon[ok] <- nm[ok]
  }

  dup <- duplicated(feature$Taxon) |
    duplicated(feature$Taxon, fromLast = TRUE)

  feature$Label <- ifelse(
    dup,
    paste0(feature$Taxon, " [", feature$TaxonID, "]"),
    feature$Taxon
  )

  feature <- feature[order(-feature$PCA_importance), ]
  top <- head(feature, Top)

  ## ---------- abundance of Top taxa ----------
  ab <- as.data.frame(rel[, top$TaxonID, drop = FALSE])
  ab$SampleID <- rownames(ab)
  ab$Group <- unname(grp[ab$SampleID])

  ab <- tidyr::pivot_longer(
    ab,
    cols = -c(SampleID, Group),
    names_to = "TaxonID",
    values_to = "Abundance"
  )

  ab <- dplyr::left_join(
    ab,
    top[, c("TaxonID", "Label", "PCA_importance")],
    by = "TaxonID"
  )

  ab_sum <- ab |>
    dplyr::group_by(TaxonID, Label, Group) |>
    dplyr::summarise(
      Mean = mean(Abundance),
      SD = stats::sd(Abundance),
      SE = SD / sqrt(dplyr::n()),
      .groups = "drop"
    )

  ## ---------- two-group statistics ----------
  gr <- unique(as.character(ab$Group))
  test <- NULL

  if (length(gr) == 2) {
    test <- dplyr::bind_rows(lapply(split(ab, ab$TaxonID), function(d) {
      x1 <- d$Abundance[d$Group == gr[1]]
      x2 <- d$Abundance[d$Group == gr[2]]

      data.frame(
        TaxonID = d$TaxonID[1],
        Group1 = gr[1],
        Group2 = gr[2],
        Mean1 = mean(x1),
        Mean2 = mean(x2),
        log2FC = log2((mean(x2) + 1e-6) / (mean(x1) + 1e-6)),
        P = tryCatch(stats::wilcox.test(x1, x2)$p.value,
                     error = function(e) NA_real_)
      )
    }))

    test$Padj <- stats::p.adjust(test$P, "BH")
    top <- dplyr::left_join(top, test, by = "TaxonID")
  }

  ## ---------- PCA variance ----------
  variance <- data.frame(
    PC = paste0("PC", seq_along(var_exp)),
    Variance = var_exp * 100,
    Cumulative = cum_exp * 100
  )

  ## ---------- sample scores ----------
  scores <- as.data.frame(pca$x)


  scores$SampleID <- rownames(scores)
  scores$Group <- unname(grp[scores$SampleID])

  p_pca <- ggplot2::ggplot(
    scores,
    ggplot2::aes(PC1, PC2, color = Group, fill = Group)
  ) +
    ggplot2::geom_point(size = 4, shape = 21) +
    ggplot2::labs(
      x = sprintf("PC1 (%.2f%%)", variance$Variance[1]),
      y = sprintf("PC2 (%.2f%%)", variance$Variance[2])
    ) +
    theme_nature()

  ## ---------- scree ----------
  nshow <- min(10, nrow(variance))

  p_scree <- ggplot2::ggplot(
    variance[1:nshow, ],
    ggplot2::aes(x = factor(PC, levels = PC[1:nshow]), y = Variance)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_point() +
    ggplot2::geom_line(ggplot2::aes(group = 1)) +
    ggplot2::labs(x = NULL, y = "Explained variance (%)") +
    ggplot2::theme_classic()

  ## ---------- importance ----------
  order_tax <- top$Label[order(top$PCA_importance)]
  top$Label <- factor(top$Label, levels = order_tax)

  p_imp <- ggplot2::ggplot(
    top,
    ggplot2::aes(Label, PCA_importance)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = Label, y = 0, yend = PCA_importance),
      linewidth = 1.5
    ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = paste0(
        "Weighted PCA importance (%)\n",
        "PC1-PC", k,
        ", cumulative = ",
        round(cum_exp[k] * 100, 1), "%"
      )
    ) +
    ggplot2::theme_classic()

  ## ---------- abundance ----------
  ab_sum$Label <- factor(ab_sum$Label, levels = order_tax)

  p_abund <- ggplot2::ggplot(
    ab_sum,
    ggplot2::aes(Label, Mean, color = Group, group = Group)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, Mean - SE), ymax = Mean + SE),
      position = ggplot2::position_dodge(0.65),
      width = 0.2
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(0.65),
      size = 3
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Relative abundance (%)") +
    ggplot2::theme_classic()

  ## ---------- combined feature figure ----------
  p_combined <- ggpubr::ggarrange(
    p_imp, p_abund,
    ncol = 2,
    widths = c(1.05, 1.35),
    align = "hv"
  )

  ## 保持旧函数前两个返回位置兼容
  list(
    p_combined,                   # [[1]]
    feature,                      # [[2]]
    pca_plot = p_pca,
    scree_plot = p_scree,
    importance_plot = p_imp,
    abundance_plot = p_abund,
    top_features = top,
    abundance_raw = ab,
    abundance_summary = ab_sum,
    group_test = test,
    pca_scores = scores,
    pca_variance = variance,
    loadings = as.data.frame(load),
    pca = pca,
    n_pc_importance = k
  )
}
